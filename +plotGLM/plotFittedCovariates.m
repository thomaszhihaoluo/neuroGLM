function plotFittedCovariates(stats,varargin)
    p=inputParser;
    p.addParameter('matchylim',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    covars = {stats.dspec.covar.label};
    matches={'left','right'};
    groups = group_covars(covars,matches);
    color.left='b';
    color.right='r';
    n_subplot_columns = ceil(sqrt(length(groups)));
    mmx=[];
    for kCov = 1:length(groups)
        count=0;
        for j=1:length(groups{kCov})     
            count=count+1;
            label = groups{kCov}{j};
            subplot(n_subplot_columns,n_subplot_columns, kCov);hold on;
            if length(groups{kCov})==1
                shadeplot2(stats.ws.(label).tr/1000, stats.ws.(label).data, sqrt(stats.wvars.(label).data));
            else
                shadeplot2(stats.ws.(label).tr/1000, stats.ws.(label).data, sqrt(stats.wvars.(label).data),'Color', color.(matches{j}));  
%                 if count>1
%                     h(count).patch.FaceAlpha=0.5;
%                     if strcmp(label,'left_clicks') || strcmp(label,'right_clicks')
%                         legend([h.mainLine],{'left clicks','right clicks'})                
%                     else
%                         legend([h.mainLine],{'left choice','right choice'})                                    
%                     end
%                 end
            end
            mmx = minmax([mmx(:);stats.ws.(label).data(:)]);            
            title_str = strrep(groups{kCov}{j},'_',' ');
            title_str = strrep(title_str,matches{j},'');
            title(title_str);
            set(gca,'xlim',minmax(stats.ws.(label).tr(:)'/1000));
            xlabel('time (s)');
        end
    end
    if params.matchylim
        matchylim(gcf,'ylim',[-1 1]*max(abs(mmx))*1.1);
    end
    set(gcf,'position',[ 1000.3          221       1101.3         1110]);
end

function groups = group_covars(covars,matches)
    count=0;
    groups={};
    for i=1:length(covars)
        matched=false(1,length(matches));        
        if ismember(covars{i},[groups{:}])
            continue
        end
        for k=1:length(matches)
            if contains(covars{i},matches{k})
                matched(k)=true;
                count=count+1;
                groups{count}=covars(i);
                other_match_patterns = setdiff(matches,matches{k});
                for j=1:length(other_match_patterns)
                    other_matches = ismember(covars,strrep(covars{i},matches{k},other_match_patterns{j}));
                    groups{count} = unique([groups{count} covars(other_matches)]);
                end
                break
            end
            if ~any(matched) && k==length(matches)
                count=count+1;
                groups{count}=covars(i);        
            end
        end     
    end 
end