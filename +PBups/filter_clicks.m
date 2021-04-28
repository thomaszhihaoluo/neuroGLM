function click_amplitudes = filter_clicks(t_clicks_s, phi, tau)
% FILTER_CLICKS Adapt or facilitate clicks
% FILTERED_CLICKS = FILTER_CLICKS(T_CLICKS, PHI, TAU) 
%
% t_clicks_s is a matrix of click times of size m X n where m is the number
% of clicks and n is the number of trials. The array is NaN-padded for
% all trials except the ones with the longest number of clicks.
% adaptation occurrs between clicks independent of which side they
% occurred on.
%
% Note that while you can have NaNs to pad the matrix, calculation
% assumes that NaNs are all at the end of each column. If this
% assumption is broken, all bups after a NaN in the same column are
% ignored. One could add a check for this, but it would be expensive.
validateattributes(phi, {'numeric'}, {'>=', 0});
validateattributes(tau, {'numeric'}, {'>=', 0});
adapted=double(~isnan(t_clicks_s));
if phi==1 || tau==0
    return
else
    ici = diff(t_clicks_s);
    expterm = exp(-ici/tau);
    for i = 2:size(t_clicks_s,1)
        %last = tau * log(1 - adapted(i-1,:)*phi);
        %check = 1 - exp((-ici(i-1,:) + last)/tau);
        adapted(i,:) = 1 + expterm(i-1,:).*(adapted(i-1,:)*phi-1);        
    end        
end
% t1 = floor(min(t_clicks_s)/params.bin_size_s)*params.bin_size_s;
% t2 = ceil(max(t_clicks_s)/params.bin_size_s)*params.bin_size_s;
% t_steps_s = t1:params.bin_size_s:t2; % row
% t_steps_s = repmat(t_steps_s, numel(t_clicks_s), 1);
% t_clicks_s = t_clicks_s(:); %column
% t_steps_s(t_steps_s<t_clicks_s)=NaN;
% if phi > 1
%     click_amplitudes = phi*nansum(exp( -(t_steps_s-t_clicks_s)./tau ));
% else
%     click_amplitudes = 1- (1-phi)*nansum(exp( -(t_steps_s-t_clicks_s)./tau ));
% end
% plot(t_steps_s, click_amplitudes)
% hold on
% for t = 1:numel(t_clicks_s)
%     plot(t_clicks_s*[1,1], [1,1], 'ko')
% end