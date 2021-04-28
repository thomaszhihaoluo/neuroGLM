function bases = makeNonlinearRaisedCos(nBases, binSize, endPoints, nlOffset, varargin)
% Make nonlinearly stretched basis consisting of raised cosines.
% Nonlinear stretching allows faster changes near the event.
%
% 	nBases: [1] - # of basis vectors
%	binSize: time bin size (separation for representing basis
%   endPoints: [2 x 1] = 2-vector containg [1st_peak  last_peak], the peak
%          (i.e. center) of the last raised cosine basis vectors
%   nlOffset: [1] offset for nonlinear stretching of x axis:  y = log(t+nlOffset)
%         (larger nlOffset -> more nearly linear stretching)
%
%  Outputs:  iht = time lattice on which basis is defined
%            ihbasis = basis itself
%            ihctrs  = centers of each basis function
%
%  Example call
%  bases = basisFactory.makeNonlinearRaisedCos(10, 1, [0 500], 2);
%
%   begin_at_0: if TRUE, the basis at endPoints(1:2) must equal 0
input_parser = inputParser;
addParameter(input_parser, 'begin_at_0', false, @(x) isscalar(x) && islogical(x))
parse(input_parser, varargin{:});

% nonlinearity for stretching x axis (and its inverse)
nlin = @(x)(log(x + 1e-20));
invnl = @(x)(exp(x) - 1e-20);

if nlOffset <= 0
    error('nlOffset must be greater than 0');
end
ff = @(x,c,dc) (cos(max(-pi, min(pi, (x-c)*pi/dc/2))) + 1)/2; % the period is 4*dc
for p=1:2
    if endPoints(1)>=0 && p==1
        mxt=endPoints(2);
        yrnge = nlin(endPoints + nlOffset);        
        yrnge(2) = (yrnge(2).*(nBases-1) + 2*yrnge(1)) ./ (nBases+1);
%             db = diff(yrnge) / (nBases-1); % spacing between raised cosine peaks
%             ctrs = yrnge(1):db:yrnge(2); % centers for basis vectors
        [db, ctrs] = calc_spacing_centers(yrnge, nBases, input_parser.Results.begin_at_0);
        iht = (0:binSize:mxt)';
        ihbasis = ff(repmat(nlin(iht + nlOffset), 1, nBases), repmat(ctrs, numel(iht), 1), db);
        ihctrs = invnl(ctrs);        
    elseif endPoints(1)<0 && endPoints(2)>0
        if p==1
            fracNeg = abs(endPoints(1))./diff(endPoints);
            nb = floor(fracNeg*(nBases+1));
            mxt=abs(endPoints(1));
            yrnge = nlin([0 abs(endPoints(1))] + nlOffset);            
            yrnge(2) = (yrnge(2).*(nb-1) + 2*yrnge(1)) ./ (nb+1);            
%             db = diff(yrnge) / (nb-1); % spacing between raised cosine peaks
%             ctrs = yrnge(1):db:yrnge(2); % centers for basis vectors
            [db, ctrs] = calc_spacing_centers(yrnge, nb, input_parser.Results.begin_at_0);
            iht1 = (0:binSize:mxt)';
            ihbasis = rot90(ff(repmat(nlin(iht1 + nlOffset), 1, nb), repmat(ctrs, numel(iht1), 1), db),2);
            ihctrs = -invnl(ctrs);             
        else
            nb1=nb;
            nb = nBases+1-nb;           
            yrnge = nlin([0 endPoints(2)] + nlOffset);
            mxt=endPoints(2);
            yrnge(2) = (yrnge(2).*(nb-1) + 2*yrnge(1)) ./ (nb+1);             
%             db = diff(yrnge) / (nb-1); % spacing between raised cosine peaks
%             ctrs = yrnge(1):db:yrnge(2); % centers for basis vectors
            [db, ctrs] = calc_spacing_centers(yrnge, nb, input_parser.Results.begin_at_0);
            iht = (0:binSize:mxt)';            
            ihbasis2=ff(repmat(nlin(iht + nlOffset), 1, nb), repmat(ctrs, numel(iht), 1), db);
            ihbasis = blkdiag(ihbasis,ihbasis2);
            ihbasis(:,nb1) = sum(ihbasis(:,[nb1 nb1+1]),2);
            ihbasis(:,nb1+1)=[];
            ihctrs = [ihctrs invnl(ctrs)];  
            iht = [-iht1;(0:binSize:mxt)'];            
        end
    elseif p==1
        error('Endpoints much be monotonically increasing and the second must be positive');
    end
end

bases.type = mfilename;
bases.param.nBases = nBases;
bases.param.binSize = binSize;
bases.param.endPoints = endPoints;
bases.param.nlOffset = nlOffset;
bases.B = ihbasis;
bases.edim = size(bases.B, 2);
bases.tr = sort(iht);
bases.centers = sort(ihctrs);

end

function [db, ctrs] = calc_spacing_centers(yrnge, nBases, begin_at_0)
    db = diff(yrnge) / (nBases-1); % spacing between raised cosine peaks
    if begin_at_0
        yrnge = yrnge + db*[2,0];
    end
    db = diff(yrnge) / (nBases-1); % spacing between raised cosine peaks
    ctrs = yrnge(1):db:yrnge(2); % centers for basis vectors
end