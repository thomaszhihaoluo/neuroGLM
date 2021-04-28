
function [h,filtered_spikes] = plotPETH(expt,variable,trialIndices,varargin)
    p=inputParser;
    p.addParameter('align_to','');
    p.addParameter('desired_delay',1000);    
    p.addParameter('bin_size_s',0.1);
    p.parse(varargin{:});
    params=p.Results;
    sd=50;
    count=0;
    trialIndices = trialIndices(:)';
    for k=trialIndices
        count=count+1;
        if ischar(variable)
            spikes = buildGLM.getBinnedSpikeTrain(expt, variable, k)*1000;
        else
            spikes = variable(buildGLM.getSpikeIndicesforTrial(expt,k))*1000;
        end
        filtered_spikes(count,:) = filterArray(full(spikes),gauss_kernel(sd*5,sd));
    end
    aligned_to = expt.param.aligned_to;
    event_time = expt.trial(1).(aligned_to); 
    if ~isempty(params.align_to)
        filtered_spikes = realign_spikes(expt,filtered_spikes,trialIndices,params.align_to,params.desired_delay);
        aligned_to = params.align_to;
        event_time = params.desired_delay;
    end
    times = ((1:expt.trial(1).duration)-event_time)/1000;
    if numel(trialIndices)>1
        min_trials =10;
        idx = sum(~isnan(filtered_spikes),1)>min_trials;
        times = times(idx);
        filtered_spikes = filtered_spikes(:,idx);
        n_bins=max(ceil(range(times)/params.bin_size_s),1);
        [a,b] = discretize(times,n_bins);
        bin_centers=mean([b(1:end-1);b(2:end)]);
        for i=1:n_bins
            filtered_binned_spikes(:,i) = nanmean(filtered_spikes(:,find(a==i)),2);
        end
        h = struct;
        h.patch = shadeplot(bin_centers,nanmean(filtered_binned_spikes),sem(filtered_binned_spikes));
        h.main_line = plot(bin_centers,nanmean(filtered_binned_spikes));
        xlabel(['Time (s) relative to ',fix_underscore(aligned_to)]);
        ylabel('Firing Rate (sp/s)');        
    else
        h=plot(times,filtered_spikes);
        xlabel(['Time (s) relative to ',fix_underscore(aligned_to)]);
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

