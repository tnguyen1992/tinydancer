function [data]=Trinh_addTriggerChannel(data,channels_to_interpolate,baby)
% extracted from Giac_EEGinterpolate
% adds a channel that didn't exist before
% created for adding Trigger Channel for MUSICOM WP4.1
% 12/10/22 Trinh Nguyen (IIT)


% microphone time-series
[event]                 = ft_read_event(baby);
event_values            = {event.value}; % Get all the event values
g                       = strcmp(event_values,'S128'); % Find S128 code
event                   = event(g); % Select only S128 rows
sample                  = [event.sample].'; % extract sample info as numeric array
triggerts               = zeros(length(data.trial{1, 1} ),1);
for j = sample'
    triggerts(j, 1) = 1 ;
end

% create the channel
for ch_to_add = 1: length(channels_to_interpolate)
    tmp_ch = channels_to_interpolate{ch_to_add};    
    if isempty(find(ismember(data.label,tmp_ch)==1))==0 % channel already exists
        display(['GIAC: channel ' tmp_ch{1,1} ' already exists and will not be created']);
    elseif isempty(find(ismember(data.label,tmp_ch)==1))==1 % channel does not exists  
        data.label{end+1,1} = tmp_ch;        
        for tr = 1: size(data.trial,2)
%             tmp_data = nan(1,size(data.trial{tr},2));
            tmp_data = triggerts';
            data.trial{1,tr}(length(data.label),:) = tmp_data;
        end
    end
end
end





