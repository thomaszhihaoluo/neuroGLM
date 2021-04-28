function [paired_kernels, paired_index] = identify_paired_kernels(kernel_names, varargin)
p=inputParser;
p.addParameter('pairing',{'left', 'right'}, @(x) validateattributes(x, {'cell'}, {'numel', 2}));
p.parse(varargin{:});
P=p.Results;

validateattributes(kernel_names, {'cell', 'string'}, {})

kernel_names = string(kernel_names);

inds_1 = find(contains(kernel_names, P.pairing{1}));
inds_2 = find(contains(kernel_names, P.pairing{2}));

paired_kernels = {};
paired_index = [];
k = 0;
for i = inds_1(:)'
for j = inds_2(:)'
    if strcmp(strrep(kernel_names{i}, P.pairing{1}, ''), ...
              strrep(kernel_names{j}, P.pairing{2}, ''))
        k = k + 1;
        paired_kernels{k,1} = strrep(kernel_names{i}, P.pairing{1}, '');
        if paired_kernels{k}(1) == '_'
            paired_kernels{k} = paired_kernels{k}(2:end);
        elseif paired_kernels{k}(end) == '_'
            paired_kernels{k} = paired_kernels{k}(1:end-1);
        end
        paired_index(k,1:2) = [i,j];
    end
end
end

