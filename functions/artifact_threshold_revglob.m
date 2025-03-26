function [ autoart ] = artifact_threshold_revglob(cfgT, data_in)

numOfTrl  = length(data_in.trialinfo);                                    % get number of trials in the data
artifact  = zeros(0, 2);                                                  % initialize artifact variable
artfctmap = cell(1, numOfTrl);                                            % preallocate artifact map

channel = ft_channelselection(cfgT.artfctdef.threshold.channel, ...
    data_in.label);

% Prune the available data to the channels of interest
for i = 1:numOfTrl
    data_in.trial{i} = data_in.trial{i}(ismember(data_in.label, channel), :);
end

% Compute mean and standard deviation across all trials for each channel
all_trials_data = cat(3, data_in.trial{:});
mean_data = mean(all_trials_data, 3);
% std_data = std(all_trials_data, 0, 3);

% Calculate global mean and standard deviation for each channel
% global_mean = mean(mean_data, 2);
global_std = std(mean_data, 0, 2);

% Check for standard deviation violations across entire trials
for i = 1:numOfTrl
    trial_data = data_in.trial{i};
    trial_std = std(trial_data,0,2);                                    % compute standard deviation for each channel over the entire trial
    
    % Check if any channel's std exceeds x SD from the global mean
    artfctmap{i} = trial_std > ( cfgT.artfctdef.threshold.stddev * global_std);            % find violations per channel
    
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
