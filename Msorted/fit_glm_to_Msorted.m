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
    p.addParameter('use_parallel',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.parse(varargin{:});
    params=p.Results;
    if params.use_parallel; n_max_workers = Inf; else n_max_workers = 1; end
    stats=[];
    tic;
    fprintf(['\n----------- ' mfilename '-----------\n...']);   
    %% make rawData and expt structures (i.e. add state times and spike times in a trial structure that neuroGLM wants)
    rawData = make_glm_trials_from_Msorted(Msorted); 
    expt=build_expt_for_pbups(rawData); 
    if isempty(params.cellno)
        params.cellno=1:length(Msorted.tt);
        if isempty(params.cellno)
            warning('No cells in Msorted!');
            return
        end
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
    parfor (c = 1:length(params.cellno), n_max_workers)
        tic
        S = struct;
        S.cellno = params.cellno(c);
        %% determine if cell is responsive enough to fit
        S.responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(S.cellno)])),expt.trial)>0)./rawData.nTrials;
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
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:rawData.nTrials);  
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
        [~, dev, stat_temp] = glmfit(dm.X, Y, 'poisson');
        fields_to_copy = fields(stat_temp);
        for f=1:length(fields_to_copy)
            S.(fields_to_copy{f}) = stat_temp.(fields_to_copy{f});
        end       
        % compute model predicted firing rates
        S.Yhat=glmval(S.beta,dm.X,'log',S.beta);
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
        %% Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
        if ~isempty(params.kfold)
            if S.badly_scaled
                fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
            else
                S.cvp = cvpartition(rawData.nTrials,'KFold',params.kfold);
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(buildGLM.addBiasColumn(dm), raw_weights , covariances);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(expt,trial_idx);        
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);    
                tic;
                Xtrain = {};
                Ytrain = {};
                Xtest = {};
                for i=params.kfold:-1:1
                    train_idx = getSpkIdxFun(S.cvp.training(i));
                    test_idx = getSpkIdxFun(S.cvp.test(i));
                    Xtrain{i} = dm.X(train_idx,:);
                    Ytrain{i} = Y(train_idx);
                    Xtest{i} = dm.X(test_idx,:);
                end
                options = statset('MaxIter',params.maxIter);
                dev = [];
                cv_stats = [];
                for i=1:params.kfold
                    [~, dev(i), cv_stats(i)] = glmfit(Xtrain{i}, Ytrain{i}, 'poisson','options',options);  
                end                  
                cvstats = rmfield(cvstats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
                for i=1:params.kfold
                    cv_stats(i).dev=dev(i);
                    cv_stats(i).Yhat=glmval(cv_stats(i).beta,Xtest{i},'log',cv_stats(i));
                    [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                end
                fprintf('Took %s.\n',timestr(toc));            
                S.cvstats=cv_stats;
                clear cv_stats
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