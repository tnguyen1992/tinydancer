% -------------------------------------------------------------------------
% FUNCTION which detects artifacts for entire trials without using a sliding window
% -------------------------------------------------------------------------
function [ autoart ] = artifact_threshold_rev(cfgT, data_in)

numOfTrl  = length(data_in.trialinfo);                                    % get number of trials in the data
artifact  = zeros(0, 2);                                                  % initialize artifact variable
artfctmap = cell(1, numOfTrl);                                            % preallocate artifact map

channel = ft_channelselection(cfgT.artfctdef.threshold.channel, ...
    data_in.label);

% Prune the available data to the channels of interest
for i = 1:numOfTrl
    data_in.trial{i} = data_in.trial{i}(ismember(data_in.label, channel), :);
end

% Check for standard deviation violations across entire trials
for i = 1:numOfTrl
    trial_data = data_in.trial{i};
    trial_std = std(trial_data, 0, 2);                                    % compute standard deviation for each channel over the entire trial
    
    % Check if any channel's std exceeds the threshold
    artfctmap{i} = trial_std < cfgT.artfctdef.threshold.stddev;           % find violations per channel
    
    if any(~artfctmap{i})                                                 % if any channel violates the threshold
        artifact = [artifact; data_in.sampleinfo(i, :)];                  % mark the entire trial as artifact
    end
end

autoart.artfctdef     = cfgT.artfctdef;                                   % generate output data structure
autoart.showcallinfo  = cfgT.showcallinfo;
autoart.artfctdef.threshold.artifact  = artifact;
autoart.artfctdef.threshold.artfctmap = artfctmap;
autoart.artfctdef.threshold.sliding   = 'no';

end
