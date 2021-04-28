function [] = plot_kernels_indiv(Stats, unit_id, varargin)
p=inputParser;
p.addParameter('ridge_parameter', 1);
p.addParameter('pairing',{'left', 'right'});
p.parse(varargin{:});
params=p.Results;

%% Plot the fits
p = params.ridge_parameter; % lamba ridge parameter
c = unit_id; % cell

kernel_names = fieldnames(Stats{1}(c).fits(p).wvars);

Gain = struct;
for i = 1:numel(Stats)
    % If Poisson, exponentiate the weights to get the multiplicative gain on the spike rate
    for fie = kernel_names(:)'; fie = fie{:};
        Gain(i,1).(fie).se   = Stats{i}(c).fits(p).wvars.(fie).data.^0.5;
        Gain(i,1).(fie).time = Stats{i}(c).fits(p).wvars.(fie).tr;
        Gain(i,1).(fie).est = Stats{i}(c).fits(p).ws.(fie).data;
        Gain(i,1).(fie).est_plus_se = Gain(i,1).(fie).est - Gain(i,1).(fie).se;
        Gain(i,1).(fie).est_minus_se = Gain(i,1).(fie).est + Gain(i,1).(fie).se;            
        if strcmp(Stats{1}.dspec.distribution, 'poisson')
            for fie2 = {'est', 'est_plus_se', 'est_minus_se', 'se'}; fie2 = fie2{:};
                Gain(i,1).(fie).(fie2) = exp(Gain(i,1).(fie).(fie2));
            end
        end
    end
end

% summarize the time duration of each kernel
end_points = struct;
for fie = kernel_names(:)'; fie = fie{:};
    idx = Stats{i}(c).dspec.idxmap.(fie);
    basis_type = Stats{i}(c).dspec.covar(idx).basis.type;
    if strcmp(basis_type, 'makeNonlinearRaisedCos')
       end_points.(fie) = Stats{i}(c).dspec.covar(idx).basis.param.endPoints;
    else
        covariate = Stats{i}(c).dspec.covar(idx);
        end_points.(fie) = covariate.offset + [0, covariate.basis.param.duration];
    end
end
%% Plot
kColor = PB_get_constant('color');
kLine = PB_get_constant('line');
[paired_kernels, paired_index] = plotGLM.identify_paired_kernels(kernel_names); % the default pairing is left, right

big_figure
n_row = 3;
n_col = 5;
p = 0;
kernel_names = fieldnames(Gain);
for k = 1:numel(kernel_names)
    kn = kernel_names{k};
    if kn == "spike_history" || any(paired_index(:) == k)
        continue
    end
    p = p + 1; subplot(n_row,n_col,p); fig_prepare_axes
    for i = 1:numel(Stats)
        shadeplot(Gain(i).(kn).time, Gain(i).(kn).est_plus_se, Gain(i).(kn).est_minus_se, 'Color', kColor.opto(i,:))
        plot(Gain(i).(kn).time, Gain(i).(kn).est, '-', 'Color', kColor.opto(i,:))
    end
    if p == 1
        title([fix_underscore(Stats{1}.dspec.expt.id) ' unit ' num2str(Stats{1}.cellno)])
    end
    quickplot('1x')
    quickplot('0y')
    xlabel(['Time from ' fix_underscore(kn) ' (' Stats{i}(c).dspec.expt.unitOfTime ')'])
    xlim(end_points.(kn))
    if strcmp(Stats{1}.dspec.distribution, 'poisson')
        ylabel('Multiplicative gain on spike rate')
    elseif strcmp(Stats{1}.dspec.distribution, 'normal')
        ylabel('Additive weight on firing rate')
    else
        ylabel('Unknown')
    end
    if contains(kn, 'click')
        set(gca, 'XScale', 'log')
    end
end
%% Plot paired kernels
for k = 1:size(paired_index,1)
    p = p + 1; subplot(n_row,n_col,p); fig_prepare_axes
    for i = 1:numel(Stats)
    for j = 1:2 % pairing
        kn = kernel_names{paired_index(k,j)};
        shadeplot(Gain(i).(kn).time, Gain(i).(kn).est_plus_se, Gain(i).(kn).est_minus_se, 'Color', kColor.opto(i,:))
        plot(Gain(i).(kn).time, Gain(i).(kn).est, kLine.default{j}, 'Color', kColor.opto(i,:))
    end
    end
    quickplot('1x')
    quickplot('0y')
    xlabel(['Time from ' fix_underscore(paired_kernels{k}) ' (' Stats{i}(c).dspec.expt.unitOfTime ')'])
    xlim(end_points.(kn))
    if strcmp(Stats{1}.dspec.distribution, 'poisson')
        ylabel('Multiplicative gain on spike rate')
    elseif strcmp(Stats{1}.dspec.distribution, 'normal')
        ylabel('Additive weight on firing rate')
    else
        ylabel('Unknown')
    end
    if contains(kn, 'click')
        set(gca, 'XScale', 'log')
    end
end