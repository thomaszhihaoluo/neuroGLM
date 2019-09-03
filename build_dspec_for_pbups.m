function dspec = build_dspec_for_pbups(expt,cellno,varargin)

    %% Build 'designSpec' which specifies how to generate the design matrix
    % Each covariate to include in the model and analysis is specified.
    dspec = buildGLM.initDesignSpec(expt);
    binfun = expt.binfun;
    left_cond = @(trial) (~trial.pokedR);
    right_cond = @(trial) (trial.pokedR);

    %% Spike history
    dspec = buildGLM.addCovariateSpiketrain(dspec, 'spike_history', ['sptrain',num2str(cellno)], 'History filter');

    %% cpoke_in
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 2500, 8, binfun);
    %bs=basisFactory.makeNonlinearRaisedCos(7,1,[-2000 1800],500);    
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in','cpoke_in', 'Acausal filter aligned to cpoke_in on left choice trials', bs,-1250);
    
    %% cpoke_in_right
    %bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 2500, 8, binfun);    
    %bs=basisFactory.makeNonlinearRaisedCos(7,1,[-2000 1800],350);    
    %dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in_right','cpoke_in', 'Acausal filter aligned to cpoke_in on right choice trials', bs,-1000,right_cond);    
    
    %% cpoke_out_left
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 5, binfun);        
    %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);    
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_left','cpoke_out', 'Acausal filter aligned to cpoke_out on left choice trials', bs,-750,left_cond);
    
    %% cpoke_out_right
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 5, binfun);            
    %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1500 1500],200);    
    dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_right','cpoke_out', 'Acausal filter aligned to cpoke_out on right choice trials', bs,-750,right_cond);    
    
    %% clicks_on_left
    bs=basisFactory.makeNonlinearRaisedCos(6,1,[0 200],10);
    dspec = buildGLM.addCovariateTiming(dspec, 'sound_on_left','clicks_on', 'Causal filter aligned to clicks_on on left choice trials', bs,0,left_cond);
    
    %% clicks_on_right
    bs=basisFactory.makeNonlinearRaisedCos(6,1,[0 200],10);
    dspec = buildGLM.addCovariateTiming(dspec, 'sound_on_right','clicks_on', 'Causal filter aligned to clicks_on on right choice trials', bs,0,right_cond);     
    
    %% left_clicks
    bs=basisFactory.makeNonlinearRaisedCos(6,1,[0 200],10);
    dspec = buildGLM.addCovariateTiming(dspec, 'left_clicks',[], 'Causal filter aligned to clicks_on on left choice trials', bs,0,left_cond);
    
    %% right_clicks
    bs=basisFactory.makeNonlinearRaisedCos(6,1,[0 200],10);
    dspec = buildGLM.addCovariateTiming(dspec, 'right_clicks',[], 'Causal filter aligned to clicks_on on right choice trials', bs,0,right_cond);  
    
    %% spoke_left_hit
    bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);                
    cond = @(trial) (~trial.pokedR & trial.ishit);   
    dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_hit','spoke', 'Acausal filter aligned to spoke on left choice hit trials', bs,-500,cond);
    
    %% spoke_right_hit
    bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);               
    cond = @(trial) (trial.pokedR & trial.ishit);        
    dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_hit','spoke', 'Acausal filter aligned to spoke on right choice hit trials', bs,-500,cond); 
    
    %% spoke_left_miss
    bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);  
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);                
    cond = @(trial) (~trial.pokedR & ~trial.ishit);    
    dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_miss','spoke', 'Acausal filter aligned to spoke on left choice miss trials', bs,-500,cond);
    
    %% spoke_right_miss
    bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);   
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, binfun);                
    cond = @(trial) (trial.pokedR & ~trial.ishit);        
    dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_miss','spoke', 'Acausal filter aligned to spoke on right choice miss trials', bs,-500,cond);     

end