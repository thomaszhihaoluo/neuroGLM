function dspec = build_dspec_for_pbups(expt,cellno)

    %% Build 'designSpec' which specifies how to generate the design matrix
    % Each covariate to include in the model and analysis is specified.
    dspec = buildGLM.initDesignSpec(expt);
    binfun = expt.binfun;
    chose_left = @(trial) (~trial.pokedR);
    chose_right = @(trial) (trial.pokedR);
    %% Spike history
    bs = basisFactory.makeNonlinearRaisedCos(6, dspec.expt.binSize, [0 250], 1);
    dspec = buildGLM.addCovariateSpiketrain(dspec, 'spike_history', ['sptrain',num2str(cellno)], 'History filter', bs);
    %% cpoke_in
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 4000, 16, binfun);
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in','cpoke_in', 'Filter aligned to cpoke_in on left choice trials', bs,-1000);
%     dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in_right','cpoke_in', 'Filter aligned to cpoke_in on left choice trials', bs,-1000, chose_right);
    %% cpoke_out
    dur_ms_pre = 1000;
    dur_ms_pos = 1000;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', dur_ms_pre+dur_ms_pos, 8, binfun);        
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_left','cpoke_out', ...
                                        'Filter aligned to cpoke_out on left choice trials', bs,-dur_ms_pre,chose_left);
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_right','cpoke_out', ...
                                        'Filter aligned to cpoke_out on right choice trials', bs,-dur_ms_pre,chose_right);    
    %% left_clicks
    bs=basisFactory.makeNonlinearRaisedCos(10,1,[0 1000],5, 'begin_at_0', true);
    dspec = buildGLM.addCovariateTiming(dspec, 'left_clicks',[], 'Causal filter aligned to clicks_on on left choice trials', bs,0,chose_left); 
    %% right_clicks
    bs=basisFactory.makeNonlinearRaisedCos(10,1,[0 1000],5, 'begin_at_0', true);
    dspec = buildGLM.addCovariateTiming(dspec, 'right_clicks',[], 'Causal filter aligned to clicks_on on right choice trials', bs,0,chose_right); 
    %% laser off
    bs=basisFactory.makeNonlinearRaisedCos(5,1,[0 1000],5, 'begin_at_0', true);
    dspec = buildGLM.addCovariateTiming(dspec, 'laser_off','laser_off', 'Causal filter aligned to beginning of the laser off ramp', bs,0);
    %% side pokes
    for make_side_poke = []
        dur_ms = 1000;
        pre_ms = 0;
        n_bases = 5;
        bs = basisFactory.makeSmoothTemporalBasis('raised cosine', dur_ms, n_bases, binfun);                
        % spoke_left_hit
        cond = @(trial) (~trial.pokedR & trial.ishit);   
        dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_hit','spoke', 'Causal filter aligned to spoke on left choice hit trials', bs,pre_ms,cond);
        % spoke_right_hit
        cond = @(trial) (trial.pokedR & trial.ishit);        
        dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_hit','spoke', 'Causal filter aligned to spoke on right choice hit trials', bs,pre_ms,cond);      
        % spoke_left_miss
        cond = @(trial) (~trial.pokedR & ~trial.ishit);    
        dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_miss','spoke', 'Causal filter aligned to spoke on left choice miss trials', bs,pre_ms,cond);
        % spoke_right_miss
        cond = @(trial) (trial.pokedR & ~trial.ishit);        
        dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_miss','spoke', 'Causal filter aligned to spoke on right choice miss trials', bs,pre_ms,cond); 
    end
    %% sound on/stereo click
    for make_sound_on = 1
        % clicks_on_left
        bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);               
        dspec = buildGLM.addCovariateTiming(dspec, 'sound_on_left','clicks_on', 'Filter aligned to clicks_on on left choice trials', bs,-800,chose_left);
        % clicks_on_right
        bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);               
        dspec = buildGLM.addCovariateTiming(dspec, 'sound_on_right','clicks_on', 'Filter aligned to clicks_on on right choice trials', bs,-800,chose_right); 
    end
end