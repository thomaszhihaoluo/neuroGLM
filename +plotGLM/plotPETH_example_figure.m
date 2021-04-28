function plotPETH_example_figure(stats,cellno, varargin)
    p=inputParser;
    p.addParameter('align_to',[]);
    p.parse(varargin{:});
    params=p.Results;
    
    subplot(2,4,[1 2 5 6]);
    
    idx = arrayfun(@(x) x.cellno == cellno, stats);
    
    poked_r = [stats(idx).dspec.expt.trial.pokedR];
    
    Colors = PB_get_constant('color');
        
    hold on;
    % predictions
    [h(1),fitted_spikes_right]=plotGLM.plotPETH(stats(idx).dspec.expt, stats(idx).fits.Yhat,find(poked_r), 'align_to',params.align_to);
    h(2) = plotGLM.plotPETH(stats(idx).dspec.expt,stats(idx).fits.Yhat,find(~poked_r), 'align_to',params.align_to);    
    %data
    [h(3),observed_spikes_right]= plotGLM.plotPETH(stats(idx).dspec.expt,['sptrain',num2str(cellno)],find(poked_r), 'align_to',params.align_to);
    h(4) = plotGLM.plotPETH(stats(idx).dspec.expt,['sptrain',num2str(cellno)],find(~poked_r), 'align_to',params.align_to);
    h(1).main_line.Color = 'r';
    h(2).main_line.Color = 'b';
    h(3).main_line.Color = 'r';
    h(4).main_line.Color = 'b';
    for i=1:4
        h(i).patch.FaceAlpha = 0;
    end
    for i=1:2
        h(i).main_line.LineStyle = '--';
    end
    correl = corr(fitted_spikes_right(:),observed_spikes_right(:));
    legend([h([3 4 1 2]).main_line],{'observed, right choice','observed, left choice','fit, right choice','fit, left choice'});
    fig_yMin0
    title('Trial average');    
    rand_trials = randi(length(stats(idx).dspec.expt.trial),4,1);
    clear h;
    subplots = [3 4 7 8];
    for i=1:4
        subplot(2,4,subplots(i));
        %pred
        h(1)=plotGLM.plotPETH(stats(idx).dspec.expt,stats(idx).fits.Yhat,rand_trials(i));    hold on;
        %data
        yyaxis right
        h(2)=plotGLM.plotPETH(stats(idx).dspec.expt,['sptrain',num2str(cellno)],rand_trials(i));
        h(1).Color=zeros(1,3);
        h(2).Color=Colors.default(2,:);
        legend(h([2 1]),{'data','fit'});
        title(['Trial ',num2str(rand_trials(i))]);
    end
    set(gcf,'position',[72.333       239.67       2374.7       751.33]);
end