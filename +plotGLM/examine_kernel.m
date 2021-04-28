function [] = examine_kernel(Stats, varargin)
p=inputParser;
p.addParameter('condition_name',{});
p.addParameter('kernel','cpoke_out');
p.addParameter('comparison',{'left', 'right'});
p.addParameter('unit_id',[]);
p.parse(varargin{:});
params=p.Results;

Kernels = plotGLM.extract_kernels(Stats);
kernel_names = fieldnames(Kernels);
n_cond = numel(Stats);
n_units = numel(Stats{1});
unit_id = params.unit_id;
if isempty(unit_id)
    unit_id = repmat("", n_units,1);
end
group_id = unique(unit_id);
n_groups = numel(group_id);
%% Extract kernel difference
k1 = contains(kernel_names, params.kernel) & contains(kernel_names, params.comparison{1});
k2 = contains(kernel_names, params.kernel) & contains(kernel_names, params.comparison{2});
if sum(k1) == 1 || sum(k2) == 1
    k1 = kernel_names{k1};
    k2 = kernel_names{k2};
    for i = 1:n_cond
        kernel_diff{i} = Kernels.(k2).ws{i} - Kernels.(k1).ws{i};
    end
    tr = Kernels.(k1).tr{i};
else % treat cells from different kernel groups as though they are the same cells
    for i = 1:n_cond
        k1_ws = [];
        k2_ws = [];
        for j =find(k1')
            k1_ws = [k1_ws; Kernels.(kernel_names{j}).ws{i}];
        end
        for j =find(k2')
            k2_ws = [k2_ws; Kernels.(kernel_names{j}).ws{i}];
        end
        kernel_diff{i} = k2_ws - k1_ws;
    end
    tr = Kernels.(kernel_names{j}).tr{i};
    unit_id = repmat(unit_id, sum(k1),1);
end
%% Take different time components
for i = 1:n_cond
    if params.kernel == "cpoke_out" || params.kernel == "sound_on"
        kernel_diff{i} = kernel_diff{i}(:, tr<=0);
    elseif contains(params.kernel, 'spoke') || contains(params.kernel, 'click')
        kernel_diff{i} = kernel_diff{i}(:, tr>0);
    end
end
%% Get either peak or sum across kernels    
for i = 1:n_cond
    if params.kernel == "cpoke_out" || params.kernel == "sound_on" || params.kernel == "spoke"
        % take peak
        [~, ind_peak] = max(abs(kernel_diff{i}), [], 2);
        ker_diff{i} = nan(numel(ind_peak),1);
        for n = 1:numel(ind_peak)
            ker_diff{i}(n,1) = kernel_diff{i}(n,ind_peak(n));
        end
    else
        ker_diff{i} = sum(kernel_diff{i}, 2); % sum across time
    end
end
    %% Plot
min_sel = 0.001;
idx_min = cellfun(@(x) abs(x) > min_sel, ker_diff, 'uni', 0);
idx_min = any(cell2mat(idx_min(:)'),2); % either condition exceedng minimum selectivity

for i = 1:n_cond
for j = i+1:n_cond
    fig_subplot(n_groups,3)
    for g = 1:n_groups
        idx = idx_min & unit_id == group_id(g);
        subplot(n_groups,2,(g-1)*2+1)
        fig_prepare_axes
        set(gca, 'DataAspectRatio', [1,1,1])
        plot(ker_diff{i}(idx), ker_diff{j}(idx), 'ko');
        title([group_id{g}, params.comparison{2} '-' params.comparison{1} ' ' fix_underscore(params.kernel) ' kernel'])
        if ~isempty(params.condition_name)
            xlabel(params.condition_name{i})
            ylabel(params.condition_name{j})
        else
            xlabel(['condition ' num2str(i)])
            ylabel(['condition ' num2str(j)])
        end
        quickplot('0xydiag')

        subplot(n_groups,2,(g-1)*2+2)    
        diff_abs_sel = abs(ker_diff{j}(idx)) - abs(ker_diff{i}(idx));
        hdl = histogram(diff_abs_sel, 20);
        legend(hdl, ['median = ' num2str(median(diff_abs_sel)), '(p = ' pval2str(signrank(diff_abs_sel)) ')'], 'Location', 'Best')
        title(sprintf('|%s| - |%s|', params.condition_name{j}, params.condition_name{i}))
        ylabel('units')
        xlabel('Kernel weight')
    end
end
end

