function stats = fit_glm_to_Msorted(Msorted,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Msorted data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('cellno',[]);
    p.addParameter('kfold',[]);
    p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_parallel',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('distribution','poisson');
    p.addParameter('expt',[]);
    p.addParameter('nTrials',[]);
    p.parse(varargin{:});
    params=p.Results;
    if params.use_parallel && length(params.cellno) > 1 
        n_max_workers = Inf; 
    else
        n_max_workers = 0; 
    end
    stats=[];
    tic;
    fprintf(['\n----------- ' mfilename '-----------\n...']);   
    %% make rawData and expt structures (i.e. add state times and spike times in a trial structure that neuroGLM wants)
    if isempty(params.expt)
        rawData = make_glm_trials_from_Msorted(Msorted, 'nTrials', nTrials); 
        expt=build_expt_for_pbups(rawData); 
        if isempty(params.cellno)
            params.cellno=1:length(Msorted.tt);
            if isempty(params.cellno)
                warning('No cells in Msorted!');
                return
            end
        end   
    else
        expt = params.expt;
    end
    %% testing save
    if params.save
        fprintf('Testing save now before spending a long time fitting... ');
        mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits_save_test');
        test=[];
        save(mat_file_name,'test','-v7.3');
        delete(mat_file_name);   
        fprintf(' Success!\n');
    end
    fprintf('\nTook %.0f seconds preparing for the glm\n', toc)
    %% loop over cells
%     parfor (c = 1:length(params.cellno), n_max_workers)
    for c = 1:length(params.cellno)
        tic
        S = struct;
        S.cellno = params.cellno(c);
        %% determine if cell is responsive enough to fit
        S.responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(S.cellno)])),expt.trial)>0)./expt.nTrials;
        if params.minResponsiveFrac>0
            if S.responsiveFrac < params.minResponsiveFrac
                fprintf('Cell %g only fired a spike on %0.2g%% of trials. Moving on without fitting.\n',S.cellno,S.responsiveFrac*100);
                for fie = fieldnames(S)'; fie = fie{:};
                    stats(c).(fie) = S.(fie);
                end
                continue
            end
        end
        %% build dspec and design matrix
        % build dspec, which states what the regressors are in your model
        % and how they are parameterized
        dspec = build_dspec_for_pbups(expt,S.cellno);
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:expt.nTrials);  
        dm = buildGLM.removeConstantCols(dm);
        S.dspec = dspec;
        Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(S.cellno)], dm.trialIndices));  
        %% determine if spike/parameter ratio is acceptable
        S.totalSpks = sum(Y);
        S.spkParamRatio = S.totalSpks ./ (size(dm.X,2)+1);   
        if params.minSpkParamRatio>0
            if S.spkParamRatio < params.minSpkParamRatio
                fprintf('Cell %g only has %g spikes to %g params to be fit. Moving on without fitting.\n',S.cellno,S.totalSpks,size(dm.X,2)+1);
                for fie = fieldnames(S)'; fie = fie{:};
                    stats(c).(fie) = S.(fie);
                end
                continue
            end
        end        
        %% Fitting UN cross-validated model
        [~, dev, stat_temp] = glmfit(dm.X, Y, params.distribution);
        toc
        fields_to_copy = fields(stat_temp);
        for f=1:length(fields_to_copy)
            S.(fields_to_copy{f}) = stat_temp.(fields_to_copy{f});
        end       
        % compute model predicted firing rates
        if strcmp(params.distribution, 'poisson')
            S.Yhat=glmval(S.beta,dm.X,'log',S.beta);
        elseif strcmp(params.distribution, 'normal')
            S.Yhat = [ones(size(dm.X,1),1), dm.X]*S.beta;
        end
        S.dev=dev;
        % reconstruct fitted kernels by weighted combination of basis functions
        [S.ws,S.wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), S.beta , S.covb);
        % determine if least-squared weights are badly scaled. If so, not much point
        % doing cross-validation.
        sqrtw=sqrt(S.wts);
        if any(sqrtw~=0 & sqrtw<(max(sqrtw)*eps('double')^(2/3)))
            S.badly_scaled=true;
        else
            S.badly_scaled=false;
        end
        S.dm = dm;
        %% Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
        if ~isempty(params.kfold)
            if S.badly_scaled
                fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
            else
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold); tic;
                cvp = cvpartition(expt.nTrials,'KFold',params.kfold);
                dm_X = dm.X; % so that d isn't broadcasted
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(buildGLM.addBiasColumn(dm), raw_weights , covariances);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(expt,trial_idx);        
                options = statset('MaxIter',params.maxIter);
                betas = S.beta;
                Xtrain = {};
                Ytrain = {};
                Xtest = {};
                parfor i=1:params.kfold                    
                    fprintf(' %i',i)
                    train_idx = getSpkIdxFun(cvp.training(i));
                    test_idx = getSpkIdxFun(cvp.test(i));
                    Xtrain{i} = dm_X(train_idx,:);
                    Ytrain{i} = Y(train_idx);
                    Xtest{i} = dm_X(test_idx,:);
                    [~, cv_dev, cv_stats] = glmfit(Xtrain{i}, Ytrain{i}, params.distribution,'options',options);
                    cv_stats.dev=cv_dev;
                    if strcmp(params.distribution, 'poisson')
                        cv_stats.Yhat=glmval(cv_stats.beta, Xtest{i},'log', cv_stats);
                    elseif strcmp(params.distribution, 'normal')
                        cv_stats.Yhat=[ones(size(Xtest{i},1),1), Xtest{i}]*cv_stats.beta;
                    end
                    [cv_stats.ws, cv_stats.wvars] = combineWeightFun(cv_stats.beta, cv_stats.covb);
                    cv_stats.mse = mean(cv_stats.resid.^2);
                    cv_stats = rmfield(cv_stats, {'resid', 'residp', 'resida' 'residd'});
                    CV_stats(i) = cv_stats;
                end
                S.cross_valid.partition = cvp;
                S.cross_valid.stats = CV_stats;
                fprintf('Took %s.\n',timestr(toc));
            end
        end
        %% Transfer information from temporary structure to the output structure
        for fie = fieldnames(S)'; fie = fie{:};
            stats(c).(fie) = S.(fie);
        end
        fprintf('\nDone with cell %i after %.0f seconds', S.cellno, toc)
    end
    %% Save
    if params.save && isfield(stats,'dev')
        mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits');
        save(mat_file_name,'stats','-v7.3');
        fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
    end
    fprintf('   took %s.\n',timestr(toc));
end