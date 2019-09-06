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
    stats=[];    
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
    fprintf('\n');
    %% testing save
    if params.save
        fprintf('Testing save now before spending a long time fitting... ');
        mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits_save_test');
        test=[];
        save(mat_file_name,'test','-v7.3');
        delete(mat_file_name);   
        fprintf(' Success!\n');
    end    
    cell_count=0;    
    for c=length(params.cellno):-1:1     % loop over cells
        cell_count=cell_count+1;
        stats(c).cellno = params.cellno(c);          
        fprintf('Cell id %g (%g of %g to fit):\n',stats(c).cellno,cell_count,length(params.cellno));
        %% determine if cell is responsive enough to fit
        stats(c).responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(stats(c).cellno)])),expt.trial)>0)./rawData.nTrials;        
        if params.minResponsiveFrac>0
            if stats(c).responsiveFrac < params.minResponsiveFrac
                fprintf('Cell %g only fired a spike on %0.2g%% of trials. Moving on without fitting.\n',stats(c).cellno,stats(c).responsiveFrac*100);
                continue
            end
        end
        %% build dspec and design matrix
        % build dspec, which states what the regressors are in your model
        % and how they are parameterized
        stats(c).dspec = build_dspec_for_pbups(expt,stats(c).cellno);    % NOTE: this should be taken out of the loop, since it's the same for all cells -AGB, 9/6/2019
        dm = buildGLM.compileSparseDesignMatrix(stats(c).dspec, 1:rawData.nTrials);  
        dm = buildGLM.removeConstantCols(dm);
        Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(stats(c).cellno)], dm.trialIndices));  
        %% determine if spike/parameter ratio is acceptable
        stats(c).totalSpks = sum(Y);
        stats(c).spkParamRatio = stats(c).totalSpks ./ (size(dm.X,2)+1);   
        if params.minSpkParamRatio>0
            if stats(c).spkParamRatio < params.minSpkParamRatio
                fprintf('Cell %g only has %g spikes to %g params to be fit. Moving on without fitting.\n',stats(c).cellno,stats(c).totalSpks,size(dm.X,2)+1);
                continue
            end
        end        
        %% Fitting UN cross-validated model
        tic;fprintf('   Fitting UN cross-validated model ... ');                    
        [~, dev, stat_temp] = glmfit(dm.X, Y, 'poisson');
        fields_to_copy = fields(stat_temp);
        for f=1:length(fields_to_copy)
            stats(c).(fields_to_copy{f}) = stat_temp.(fields_to_copy{f});
        end       
        % compute model predicted firing rates
        stats(c).Yhat=glmval(stats(c).beta,dm.X,'log',stats(c).beta);
        fprintf('took %s.\n',timestr(toc));
        stats(c).dev=dev;
        % reconstruct fitted kernels by weighted combination of basis functions
        [stats(c).ws,stats(c).wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), stats(c).beta , stats(c).covb);
        % determine if least-squared weights are badly scaled. If so, not much point
        % doing cross-validation.
        sqrtw=sqrt(stats(c).wts);
        if any(sqrtw~=0 & sqrtw<(max(sqrtw)*eps('double')^(2/3)))
            stats(c).badly_scaled=true;
        else
            stats(c).badly_scaled=false;
        end
        %% Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
        if ~isempty(params.kfold)
            if stats(c).badly_scaled
                fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
            else
                stats(c).cvp = cvpartition(rawData.nTrials,'KFold',params.kfold);
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(buildGLM.addBiasColumn(dm), raw_weights , covariances);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(expt,trial_idx);        
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);    
                tic;
                for i=params.kfold:-1:1
                    train_idx = getSpkIdxFun(stats(c).cvp.training(i));
                    test_idx = getSpkIdxFun(stats(c).cvp.test(i));
                    Xtrain{i} = dm.X(train_idx,:);
                    Ytrain{i} = Y(train_idx);
                    Xtest{i} = dm.X(test_idx,:);
                end
                options = statset('MaxIter',params.maxIter);
                if params.use_parallel
                    parfor i=1:params.kfold
                        [~, dev(i), cv_stats(i)] = glmfit(Xtrain{i}, Ytrain{i}, 'poisson','options',options);  
                    end    
                else
                    for i=1:params.kfold
                        [~, dev(i), cv_stats(i)] = glmfit(Xtrain{i}, Ytrain{i}, 'poisson','options',options);  
                    end                  
                end
                cvstats = rmfield(cvstats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
                for i=1:params.kfold
                    cv_stats(i).dev=dev(i);
                    cv_stats(i).Yhat=glmval(cv_stats(i).beta,Xtest{i},'log',cv_stats(i));
                    [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                end
                fprintf('Took %s.\n',timestr(toc));            
                stats(c).cvstats=cv_stats;
                clear cv_stats
            end
        end
    end
    if params.save && isfield(stats,'dev')
        mat_file_name = strrep(Msorted.mat_file_name,'Msorted','glmfits');
        save(mat_file_name,'stats','-v7.3');
        fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
    end
end