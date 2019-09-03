function [h,filtered_spikes] = plotPETH(expt,variable,trialIndices,varargin)
    p=inputParser;
    p.addParameter('align_to','');
    p.addParameter('desired_delay',500);    
    p.addParameter('bin_size',10);
    p.parse(varargin{:});
    params=p.Results;
    sd=50;
    count=0;
    for k=trialIndices
        count=count+1;
        if ischar(variable)
            spikes = buildGLM.getBinnedSpikeTrain(expt, variable, k)*1000;
        else
            spikes = variable(buildGLM.getSpikeIndicesforTrial(expt,k))*1000;
        end
        filtered_spikes(count,:) = filterArray(full(spikes),my_gauss_kernel(sd*5,sd));
    end
    aligned_to = expt.param.aligned_to;
    event_time = expt.trial(1).(aligned_to); 
    if ~isempty(params.align_to)
        filtered_spikes = realign_spikes(expt,filtered_spikes,trialIndices,params.align_to,params.desired_delay);
        aligned_to = params.align_to;
        event_time = params.desired_delay;
    end
    min_trials =10;
    times = ((1:expt.trial(1).duration)-event_time)/1000;
    idx = sum(~isnan(filtered_spikes))>min_trials;
    if numel(trialIndices)>1


times = times(idx);
filtered_spikes = filtered_spikes(:,idx);
n_bins=100;
[a,b] = discretize(times,n_bins);
bin_centers=mean([b(1:end-1);b(2:end)]);
for i=1:n_bins
    filtered_binned_spikes(:,i) = nanmean(filtered_spikes(:,find(a==i)),2);
end


        
        
        h=shadedErrorBar(bin_centers,filtered_binned_spikes,{@nanmean,@SE});
        xlabel(['Time (s) relative to ',aligned_to]);
        ylabel('Firing Rate (sp/s)');        
    else
        h=plot(times(idx),filtered_spikes(idx));
        xlabel(['Time (s) relative to ',aligned_to]);
        ylabel('Firing Rate (sp/s)');
    end
end

function realigned_spikes = realign_spikes(expt,spikes,trialIndices,align_to,desired_delay)
    realigned_spikes = NaN(size(spikes));
    for i=1:length(trialIndices)
        event_time = expt.trial(trialIndices(i)).(align_to);
        start = event_time - desired_delay;
        start = max(1,start);
        finish = max(max(expt.trial(trialIndices(i)).left_clicks),max(expt.trial(trialIndices(i)).right_clicks));
        idx = start:finish;
        realigned_spikes(i,1:length(idx)) = spikes(i,idx);
    end
end

