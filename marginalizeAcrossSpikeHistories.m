function Y_hat = marginalizeAcrossSpikeHistories(dm,beta,covIdx)
    % covidx is the output of strcmp({dm.dspec.covar.label},'spike_history') if
    % your spike history filter is called as such
    Y_hat = zeros(size(dm.X,1),1);
    subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dm.dspec);
    covar = dm.dspec.covar(covIdx);
    basisSize = size(stats(1).dspec.covar(1).basis.B);
    sidx = subIdxs{covIdx}; 
    for t=dm.trialIndices
        spkIndices = buildGLM.getSpikeIndicesforTrial(dm.dspec.expt,t);
        for i=1:length(spkIndices)
            % this requires two steps EACH TIME POINT.        
            %% 1: equivalent of compileSparseDesignMatrix at time step t using the model-fit spike probability from time 1:t-1 as the spike history        
            if isfield(covar, 'basis') && ~isempty(covar.basis)          
%                 tmp = basisFactory.convBasis(Y_hat(spkIndices(1:i)), covar.basis, covar.offset);
%                 dm.X(spkIndices(i), sidx) = tmp(end,:);                        
                dm.X(spkIndices(i),sidx) = Y_hat(spkIndices(max(i-basisSize(1)+1,1):i))*covar.basis(end-i+1:end,:);     
            else
                %dm.X(spkIndices(i), sidx) = Y_hat_trial(i); %% analogous to
                %equivalent line in compileSparseDesignMatrix
                error('not sure how to handle this case.')
            end  
            %% 2: compute model spike probability at time t (simple glmval)        
            Y_hat(spkIndices(i))=glmval(beta,dm.X(spkIndices(i),:),'log',beta); % linkfun(x*beta); doublecheck what glmval is doing but it is something simple like this
        end
    end
end
