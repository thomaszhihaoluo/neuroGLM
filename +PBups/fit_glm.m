function stats = fit_glm(Msorted,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Msorted data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('cellno',[]);
    p.addParameter('kfold',[]);
%     p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_parallel',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('distribution','poisson');
    p.addParameter('expt',[]);
    p.addParameter('lambda_ridge',0, @(x)validateattributes(x,{'numeric'},{'vector', 'nonnegative'}));
    p.addParameter('vl_trials',[],@(x)validateattributes(x,{'logical'},{'vector'}));
    p.parse(varargin{:});
    params=p.Results;
    if params.use_parallel && length(params.cellno) > 1 
        n_max_workers = Inf; 
    else
        n_max_workers = 0; 
    end
    
    t0 = tic;
    fprintf(['\n----------- ' mfilename '-----------\n...']);   
    %% make rawData and expt structures 
    % (i.e. add state times and spike times in a trial structure that neuroGLM wants)
    if isempty(params.expt)
        rawData = PBups.make_glm_trials_from_Msorted(Msorted, 'vl_trials', params.vl_trials); 
        expt=PBups.build_expt_for_pbups(rawData); 
        if isempty(params.cellno)
            params.cellno=1:numel(Msorted.raw_spike_time_s);
            if isempty(params.cellno)
                warning('No cells in Msorted!');
                return
            end
        end   
    else
        expt = params.expt;
    end
    fprintf('\nTook %.0f seconds to prepare "rawData" and "expt" for the glm\n', toc(t0))
    %% testing save
%     if params.save
%         fprintf('Testing save now before spending a long time fitting... ');
%         mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits_save_test');
%         test=[];
%         save(mat_file_name,'test','-v7.3');
%         delete(mat_file_name);   
%         fprintf(' Success!\n');
%     end
    %% loop over cells
    parfor (c = 1:length(params.cellno), n_max_workers)
%     for c = 1:length(params.cellno)
        tic
        S = struct;
        S.cellno = params.cellno(c);
        %% determine if cell is responsive enough to fit
        S.responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(S.cellno)])),expt.trial)>0)./expt.nTrials;
        if params.minResponsiveFrac>0 && S.responsiveFrac < params.minResponsiveFrac
            fprintf('\nCell %g only fired a spike on %0.2g%% of trials. Moving on without fitting.\n',S.cellno,S.responsiveFrac*100);
            continue
        end
        %% build dspec and design matrix
        % build dspec, which states what the regressors are in your model
        % and how they are parameterized
        dspec = PBups.build_dspec_for_pbups(expt,S.cellno);
        dspec.distribution = params.distribution;
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:expt.nTrials);  
        dm = buildGLM.removeConstantCols(dm);
        S.dspec = dspec;
        Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(S.cellno)], dm.trialIndices));
        S.Y = Y;
        %% Prepare ridge regression
        Imat = eye(size(dm.X,2)+1); % identity matrix of size of filter + const
        Imat(1,1) = 0; % don't apply penalty to constant coeff
        lambda_ridge = params.lambda_ridge(:);
        n_lambda_ridge = numel(params.lambda_ridge);
        %% determine if spike/parameter ratio is acceptable
        S.totalSpks = sum(Y);
        S.spkParamRatio = S.totalSpks ./ (size(dm.X,2)+1);   
        if params.minSpkParamRatio>0 && S.spkParamRatio < params.minSpkParamRatio
            fprintf('\nCell %g only has %g spikes to %g params to be fit. Moving on without fitting.\n',S.cellno,S.totalSpks,size(dm.X,2)+1);
            continue
        end
        %% Loop over parameters for fitting
        S.lambda_ridge = params.lambda_ridge;
        for p = 1:n_lambda_ridge
            %% Fitting UN cross-validated model  
            % compute model predicted firing rates
            if strcmp(params.distribution, 'poisson')
                [~, ~, stat_temp] = glmfit(dm.X, Y, params.distribution);
                betas = stat_temp.beta;
                Yhat=glmval(betas,dm.X,'log',betas);
                residuals = stat_temp.resid;
                covb = stat_temp.covb;
            elseif strcmp(params.distribution, 'normal')
                X1 = [ones(size(dm.X,1),1), dm.X]; % design matrix with the column of 1's
                betas = ((X1'*X1)+Imat*params.lambda_ridge(p))\(X1'*Y); % ordinary least square
                Yhat = X1*betas; % predicted spike train
                residuals = Y-Yhat;
                covb = (residuals'*residuals)/(size(X1,1) - size(X1,2))*inv(X1'*X1);
            end
            % reconstruct fitted kernels by weighted combination of basis functions
            [ws,wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), betas , covb);
            S.fits(p,1).betas = betas;
            S.fits(p,1).Yhat = Yhat;
            S.fits(p,1).residuals = residuals;
            S.fits(p,1).covb = covb;
            S.fits(p,1).ws = ws;
            S.fits(p,1).wvars = wvars;
            %% Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
            if ~isempty(params.kfold) && strcmp(params.distribution, 'normal')
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);
                cvp = cvpartition(expt.nTrials,'KFold',params.kfold);
                dm_X = dm.X; % so that d isn't broadcasted
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(buildGLM.addBiasColumn(dm), raw_weights , covariances);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(expt,trial_idx);        
                options = statset('MaxIter',params.maxIter);
                cv_mse = nan(params.kfold,1);
                cv_mse1 = nan(params.kfold,1); % MSE with only a constant term
                for i=1:params.kfold                    
                    fprintf(' %i',i)
                    train_idx = getSpkIdxFun(cvp.training(i));
                    test_idx = getSpkIdxFun(cvp.test(i));
                    Xtrain = dm_X(train_idx,:);
                    Ytrain = Y(train_idx);
                    Xtest = dm_X(test_idx,:);
                    Ytest = Y(test_idx);
                    if strcmp(params.distribution, 'normal')
                        X1_train = [ones(size(Xtrain,1),1), Xtrain]; % design matrix with the column of 1's
                        cv_beta = ((X1_train'*X1_train)+Imat*params.lambda_ridge(p))\(X1_train'*Ytrain); % ordinary least square
                        X1_test = [ones(size(Xtest,1),1), Xtest];
                        cv_Yhat = X1_test*cv_beta; % predicted spike train
                        cv_resid = Ytest - cv_Yhat;
                        cv_mse(i) = mean(cv_resid.^2);
                        cv_mse1(i) = mean((Ytest - mean(Ytrain)).^2);
                    end
                end
                S.CV_mse{p,1} = cv_mse;
                S.CV_mse1{p,1} = cv_mse1;
                S.cv_mse_med = cellfun(@nanmedian, CV_mse);
                S.cv_mse1_med = cellfun(@nanmedian, CV_mse1);
            end
        end
        %% Transfer information from temporary structure to the output structure
        for fie = fieldnames(S)'; fie = fie{:};
            stats(c).(fie) = S.(fie);
        end
        fprintf('\nDone with cell %i after %.0f seconds', S.cellno, toc)
    end
    %% Save
%     if params.save %&& isfield(stats,'dev')
%         mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits');
%         save(mat_file_name,'stats','-v7.3');
%         fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
%     end
%     fprintf('   took %s.\n',timestr(toc(t0)));
end