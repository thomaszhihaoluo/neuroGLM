function [] = examine_cpoke_out(Stats, varargin)
p=inputParser;
p.addParameter('condition_name',{});
p.parse(varargin{:});
params=p.Results;

Kernels = plotGLM.extract_kernels(Stats);
kernel_names = fieldnames(Kernels);
n_cond = numel(Stats);
for i = 1:n_cond
    if any(kernel_names == "cpoke_out_left") && ...
       any(kernel_names == "cpoke_out_right")
            R_minus_L = Kernels.cpoke_out_right.ws{i} - Kernels.cpoke_out_left.ws{i};
            R_minus_L = R_minus_L(:, Kernels.cpoke_out_left.tr{i}<=0);
            [~, ind_peak] = max(abs(R_minus_L), [], 2);
            peak_R_minus_L{i} = nan(numel(ind_peak),1);
            for n = 1:numel(ind_peak)
                peak_R_minus_L{i}(n,1) = R_minus_L(n,ind_peak(n));
            end
    end
end

min_sel = 0.001;
idx = cellfun(@(x) abs(x) > min_sel, peak_R_minus_L, 'uni', 0);
idx = any(cell2mat(idx(:)'),2); % either condition exceedng minimum selectivity

for i = 1:n_cond
for j = i+1:n_cond
    fig_subplot(1,3)
    subplot(1,2,1)
    fig_prepare_axes
    set(gca, 'DataAspectRatio', [1,1,1])
    plot(peak_R_minus_L{i}(idx), peak_R_minus_L{j}(idx), 'ko');
    title('peak R-L premovement kernel')
    if ~isempty(params.condition_name)
        xlabel(params.condition_name{i})
        ylabel(params.condition_name{j})
    else
        xlabel(['condition ' num2str(i)])
        ylabel(['condition ' num2str(j)])
    end
    quickplot('0xydiag')
    
    subplot(1,2,2)    
    diff_abs_sel = abs(peak_R_minus_L{j}(idx)) - abs(peak_R_minus_L{i}(idx));
    hdl = histogram(diff_abs_sel, 20);
    legend(hdl, ['median = ' num2str(median(diff_abs_sel)), '(p = ' pval2str(signrank(diff_abs_sel)) ')'], 'Location', 'Best')
    title(sprintf('|%s| - |%s|', params.condition_name{j}, params.condition_name{i}))
end
end

