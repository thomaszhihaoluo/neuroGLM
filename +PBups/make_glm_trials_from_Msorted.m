function rawData = make_glm_trials_from_Msorted(Msorted,varargin)
    % Makes the rawData trial structure used by neuroGLM from an Msorted
    % file.
    % The final structure has fields giving the time of various events on
    % each trial and the spike times for each unit (those fields are called
    % "sptrain1",'sptrain2", etc
    %units in s
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addParameter('samplingFreq',1e3,@(x)validateattributes(x,{'numeric'},{'scalar'}));
    p.addParameter('removeViolations',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    % time relative to reference event over which you include spikes (make sure the window 
    % over which you include spiking data in your input structure is the same or smaller)    
    p.addParameter('spikeWindowS',[-2 4],@(x)validateattributes(x,{'numeric'},{'numel',2})); 
    p.addParameter('vl_trials',[],@(x)validateattributes(x,{'logical'},{'vector'}));
    
    p.parse(varargin{:});
    params = p.Results;
    fields_to_copy = {'rat','sess_date','sess_time','mat_file_name','sessid'};
    for f=1:length(fields_to_copy)
        rawData.param.(fields_to_copy{f}) = Msorted.(fields_to_copy{f});
    end
    events = {'cpoke_in','cpoke_out','clicks_on','spoke','left_clicks','right_clicks', 'laser_on', 'laser_off'};
    duration = round(diff(params.spikeWindowS)*params.samplingFreq);
    %% preallocate structure in memory
    rawData.trial = struct();
    rawData.trial(Msorted.nTrials).duration = duration; % preallocate
    n_units = numel(Msorted.raw_spike_time_s);
    
    for t = 1:Msorted.nTrials
        trial_start_time = Msorted.Trials.stateTimes.(params.ref_event)(t) + params.spikeWindowS(1);
        rawData.trial(t).duration = duration;
        rawData.trial(t).ishit = Msorted.Trials.is_hit(t);
        for e=1:length(events)
            if strcmp(events{e},'left_clicks')
               rawData.trial(t).left_clicks  = Msorted.Trials.stateTimes.left_clicks(Msorted.Trials.stateTimes.left_click_trial==t) - trial_start_time;
            elseif strcmp(events{e},'right_clicks')
               rawData.trial(t).right_clicks  = Msorted.Trials.stateTimes.right_clicks(Msorted.Trials.stateTimes.right_click_trial==t) - trial_start_time;
            else
                rawData.trial(t).(events{e}) = Msorted.Trials.stateTimes.(events{e})(t) - trial_start_time;                    
            end
            if isnan(rawData.trial(t).(events{e}))
                rawData.trial(t).(events{e}) = [];
            else
                rawData.trial(t).(events{e}) = round(rawData.trial(t).(events{e})*params.samplingFreq);
            end
        end
        rawData.trial(t).pokedR = Msorted.Trials.pokedR(t);
        for c = 1:n_units
            spike_time_s = Msorted.spike_time_s.(params.ref_event){c}{t};
            spike_time_s(spike_time_s < params.spikeWindowS(1) | spike_time_s > params.spikeWindowS(2)) = []; % i
            rawData.trial(t).(['sptrain',num2str(c)]) = round( (spike_time_s - params.spikeWindowS(1) ) * params.samplingFreq);
        end
    end
    if ~isempty(p.Results.vl_trials)
        exclude_trials= ~p.Results.vl_trials;
    end
    if params.removeViolations
        exclude_trials = exclude_trials | Msorted.Trials.violated;
    end
    rawData.nTrials = sum(~exclude_trials);
    rawData.trial = rawData.trial(~exclude_trials);    
    rawData.param.samplingFreq = params.samplingFreq; % Hz
    rawData.param.aligned_to = params.ref_event;
    rawData.param.ncells = n_units;
end