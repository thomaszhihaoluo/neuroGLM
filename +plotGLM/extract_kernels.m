function Kernels = extract_kernels(Stats, varargin)
p = inputParser;
p.addParameter('ridge_parameter', 1);
p.addParameter('vi_unit',[]);
p.parse(varargin{:});
params=p.Results;
p = params.ridge_parameter; % simplify

n_units = numel(Stats{1});
if ~isempty(params.vi_unit)
    vl_units = false(n_units,1);
    vl_units(params.vi_unit) = true;
else
    vl_units = true(n_units,1);
end

if ~iscell(Stats), Stats = {Stats}; end
Stats = cellfun(@(x) x(:), Stats, 'uni',0);
for i = 1:size(Stats,1)
for j = 1:size(Stats,2)
    kernel_names = fieldnames(Stats{i,j}(1).fits(p).ws);
    for k = 1:numel(kernel_names)
        % time vector
        tr = cell2mat(arrayfun(@(x) x.fits(p).ws.(kernel_names{k}).tr(:)', Stats{i,j}(vl_units), 'uni', 0));
        tr = unique(tr, 'rows');
        if size(tr,1) > 1; error('non-unique time vector'); end
        Kernels.(kernel_names{k}).tr{i,j} = tr;
        % weights
        Kernels.(kernel_names{k}).ws{i,j} = cell2mat(arrayfun(@(x) x.fits(p).ws.(kernel_names{k}).data(:)', Stats{i,j}(vl_units), 'uni', 0));
        Kernels.(kernel_names{k}).wvars{i,j} = cell2mat(arrayfun(@(x) x.fits(p).wvars.(kernel_names{k}).data(:)', Stats{i,j}(vl_units), 'uni', 0));
        Kernels.(kernel_names{k}).se{i,j} = Kernels.(kernel_names{k}).wvars{i,j}.^0.5;
    end
end
end
