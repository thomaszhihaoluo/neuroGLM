function rawData = make_glm_trials_from_Msorted(Msorted,varargin)
    %units in s
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addParameter('samplingFreq',1e3,@(x)validateattributes(x,{'numeric'},{'scalar'}));
    p.addParameter('removeViolations',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('removeStimTrials',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.parse(varargin{:});
    params = p.Results;
    PB_set_constants;
    fields_to_copy = {'rat','sess_date','sess_time','mat_file_name','sessid'};
    for f=1:length(fields_to_copy)
        rawData.param.(fields_to_copy{f}) = Msorted.(fields_to_copy{f});
    end
    other_events = {'cpoke_in','cpoke_out','clicks_on','spoke','left_clicks','right_clicks'};
    duration = round(diff(kSpikeWindowS.(params.ref_event))*params.samplingFreq);
        
    %% preallocate structure in memory
    rawData.trial = struct();
    rawData.trial(Msorted.nTrials).duration = duration; % preallocate

    for t = 1:Msorted.nTrials
        trial_start_time = Msorted.Trials.stateTimes.(params.ref_event)(t) + kSpikeWindowS.(params.ref_event)(1);
        rawData.trial(t).duration = duration;
        rawData.trial(t).ishit = Msorted.Trials.is_hit(t);
        for e=1:length(other_events)
            if strcmp(other_events{e},'left_clicks')
               rawData.trial(t).left_clicks  = Msorted.Trials.stateTimes.left_clicks(Msorted.Trials.stateTimes.left_click_trial==t) - trial_start_time;                   ;
            elseif strcmp(other_events{e},'right_clicks')
               rawData.trial(t).right_clicks  = Msorted.Trials.stateTimes.right_clicks(Msorted.Trials.stateTimes.right_click_trial==t) - trial_start_time;
            else
                rawData.trial(t).(other_events{e}) = Msorted.Trials.stateTimes.(other_events{e})(t) - trial_start_time;                    
            end
            if isnan(rawData.trial(t).(other_events{e}))
                rawData.trial(t).(other_events{e}) = [];
            else
                rawData.trial(t).(other_events{e}) = round(rawData.trial(t).(other_events{e})*params.samplingFreq);
            end
        end
        rawData.trial(t).pokedR = Msorted.Trials.pokedR(t);
        for c = 1:length(Msorted.tt)
            rawData.trial(t).(['sptrain',num2str(c)]) = round( (Msorted.spike_time_s.(params.ref_event){c}{t} - kSpikeWindowS.(params.ref_event)(1) ) * params.samplingFreq);
        end
    end
    if params.removeViolations
        exclude_trials = Msorted.Trials.violated;
    end
    if params.removeStimTrials
        exclude_trials = exclude_trials | Msorted.Trials.laser.isOn;
    end
    rawData.nTrials = sum(~exclude_trials);
    rawData.trial = rawData.trial(~exclude_trials);    
    rawData.param.samplingFreq = params.samplingFreq; % Hz
    rawData.param.aligned_to = params.ref_event;
    rawData.param.ncells = length(Msorted.tt);
end