function [data]=Trinh_addTriggerChannelBack(data,channels_to_interpolate,data_trigger)
% extracted from Giac_EEGinterpolate
% adds trigger channel back into the data
% created for adding Trigger Channel for MUSICOM WP4.1
% 12/10/22 Trinh Nguyen (IIT)

for ch_to_add = 1: length(channels_to_interpolate)
    tmp_ch = channels_to_interpolate{ch_to_add};    
    if isempty(find(ismember(data.label,tmp_ch)==1))==0 % channel already exists
        display(['GIAC: channel ' tmp_ch{1,1} ' already exists and will not be created']);
    elseif isempty(find(ismember(data.label,tmp_ch)==1))==1 % channel does not exists  
        data.label{end+1,1} = tmp_ch;        
        for tr = 1: size(data.trial,2)
%             tmp_data = nan(1,size(data.trial{tr},2));
            tmp_data = data_trigger.trial{1,tr};
            data.trial{1,tr}(length(data.label),:) = tmp_data;
        end
    end
end