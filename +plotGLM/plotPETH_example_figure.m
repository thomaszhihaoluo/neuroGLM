function plotPETH_example_figure(stats,cellno)
    subplot(2,4,[1 2 5 6]);
    poked_r = [stats.dspec.expt.trial.pokedR];
    signal = cellfun(@length,{stats.dspec.expt.trial.right_clicks}) ./ ...
        cellfun(@length,{stats.dspec.expt.trial.left_clicks});
    signal(signal<1) = -1./signal(signal<1);
    non_acc_trials = isnan(signal);
    cutoff = median(abs(signal(~non_acc_trials)));
    is_hit = stats.dspec.expt.trial.ishit;
    colors = rainbowColors(4);
    for i=1:4
        if i==1
            these_trials = is_hit & signal>cutoff & ~non_acc_trials;
        elseif i==2
            these_trials = is_hit & signal<cutoff & ~non_acc_trials & signal>0;
        elseif i==3
            these_trials = is_hit & abs(signal)<cutoff & ~non_acc_trials & signal<0;            
        else
            these_trials = is_hit & signal<-cutoff & ~non_acc_trials;            
        end
        h(i)=plotGLM.plotPETH(stats.dspec.expt,'sptrain12',find(these_trials),'align_to','clicks_on');hold on;            
        h(i).mainLine.Color = colors(i,:);
        h(i).patch.FaceColor = colors(i,:);
    end        
    
%     % predictions
%     [h(1),fitted_spikes_right]=plotGLM.plotPETH(stats.dspec.expt,stats.Yhat,find(poked_r),'r');hold on;
%     h(2) = plotGLM.plotPETH(stats.dspec.expt,stats.Yhat,find(~poked_r),'b');    
%     %data
%     [h(3),observed_spikes_right]= plotGLM.plotPETH(stats.dspec.expt,['sptrain',num2str(cellno)],find(poked_r),'r');
%     h(4) = plotGLM.plotPETH(stats.dspec.expt,['sptrain',num2str(cellno)],find(~poked_r),'b');
%     correl = corr(fitted_spikes_right(:),observed_spikes_right(:));
    for i=1:4
        set(h(i).patch,'facealpha',0);
    end
    for i=1:2
        set(h(i).mainLine,'linestyle',':');
    end
    legend([h([3 4 1 2]).mainLine],{'observed, right choice','observed, left choice','fit, right choice','fit, left choice'});
    title('Trial average');    
    rand_trials = randi(length(stats.dspec.expt.trial),4,1);
    clear h;
    subplots = [3 4 7 8];
    for i=1:4
        subplot(2,4,subplots(i));
        %pred
        h(1)=plotGLM.plotPETH(stats.dspec.expt,stats.Yhat,rand_trials(i));    hold on;
        %data
        h(2)=plotGLM.plotPETH(stats.dspec.expt,['sptrain',num2str(cellno)],rand_trials(i));
        h(1).Color=[0.5 0.5 0.5];
        h(2).Color=[0 0 0];
        legend(h([2 1]),{'data','fit'});
        title(['Trial ',num2str(rand_trials(i))]);
    end
    set(gcf,'position',[72.333       239.67       2374.7       751.33]);
end