% -------------------------------------------------------------------------
% SUBFUNCTION which detects artifacts by using a sliding window
% -------------------------------------------------------------------------
function [ autoart ] = artifact_sliding_threshold(cfgT, data_in)

  numOfTrl  = length(data_in.trialinfo);                                    % get number of trials in the data
  winsize   = cfgT.artfctdef.threshold.winsize * data_in.fsample / 1000;    % convert window size from milliseconds to number of samples
  artifact  = zeros(0,2);                                                   % initialize artifact variable
  artfctmap{1,numOfTrl} = [];

  channel = ft_channelselection(cfgT.artfctdef.threshold.channel, ...
              data_in.label);

  for i = 1:1:numOfTrl
    data_in.trial{i} = data_in.trial{i}(ismember(data_in.label, ...         % prune the available data to the channels of interest
                        channel) ,:);
  end

  if isfield(cfgT.artfctdef.threshold, 'range')                             % check for range violations
    for i=1:1:numOfTrl
      tmpmin = movmin(data_in.trial{i}, winsize, 2);                        % get all minimum values
      if mod(winsize, 2)                                                    % remove useless results from the edges
        tmpmin = tmpmin(:, (winsize/2 + 1):(end-winsize/2));
      else
        tmpmin = tmpmin(:, (winsize/2 + 1):(end-winsize/2 + 1));
      end

      tmpmax = movmax(data_in.trial{i}, winsize, 2);                        % get all maximum values
      if mod(winsize, 2)                                                    % remove useless results from the edges
        tmpmax = tmpmax(:, (winsize/2 + 1):(end-winsize/2));
      else
        tmpmax = tmpmax(:, (winsize/2 + 1):(end-winsize/2 + 1));
      end
      tmp = abs(tmpmin - tmpmax);                                           % estimate a moving maximum difference

      artfctmap{i} = tmp > cfgT.artfctdef.threshold.range;                  % find all violations
      [channum, begnum] = find(artfctmap{i});                               % estimate pairs of channel numbers and begin numbers for each violation
      artfctmap{i} = [artfctmap{i} false(length(channel), winsize - 1)];    % extend artfctmap to trial size
      endnum = begnum + winsize - 1;                                        % estimate end numbers for each violation
      for j=1:1:length(channum)
        artfctmap{i}(channum(j), begnum(j):endnum(j)) = true;               % extend the violations in the map to the window size
      end
      if ~isempty(begnum)
        begnum = unique(begnum);                                            % select all unique violations
        begnum = begnum + data_in.sampleinfo(i,1) - 1;                      % convert relative sample number into an absolute one
        begnum(:,2) = begnum(:,1) + winsize - 1;
        artifact = [artifact; begnum];                                      %#ok<AGROW> add results to the artifacts matrix
      end
    end
  elseif isfield(cfgT.artfctdef.threshold, 'stddev')                        % check for standard deviation violations
    for i=1:1:numOfTrl
      tmp = movstd(data_in.trial{i}, winsize, 0, 2);                        % estimate a moving standard deviation
      if mod(winsize, 2)                                                    % remove useless results from the edges
        tmp = tmp(:, (winsize/2 + 1):(end-winsize/2));
      else
        tmp = tmp(:, (winsize/2 + 1):(end-winsize/2 + 1));
      end

      artfctmap{i} = tmp > cfgT.artfctdef.threshold.stddev;                 % find all violations
      [channum, begnum] = find(artfctmap{i});                               % estimate pairs of channel numbers and begin numbers for each violation
      artfctmap{i} = [artfctmap{i} false(length(channel), winsize - 1)];    % extend artfctmap to trial size
      endnum = begnum + winsize - 1;                                        % estimate end numbers for each violation
      for j=1:1:length(channum)
        artfctmap{i}(channum(j), begnum(j):endnum(j)) = true;               % extend the violations in the map to the window size
      end
      if ~isempty(begnum)
        begnum = unique(begnum);                                            % select all unique violations
        begnum = begnum + data_in.sampleinfo(i,1) - 1;                      % convert relative sample number into an absolute one
        begnum(:,2) = begnum(:,1) + winsize - 1;
        artifact = [artifact; begnum];                                      %#ok<AGROW> add results to the artifacts matrix
      end
    end
  elseif isfield(cfgT.artfctdef.threshold, 'mad')                           % check for median absolute deviation violations
    data_continuous = cat(2, data_in.trial{:});                             % concatenate all trials
    tmpmad = mad(data_continuous, 1, 2);                                    % estimate the median absolute deviation of the whole data
    tmpmedian = median(data_continuous, 2);                                 % estimate the median of the data

    for i=1:1:numOfTrl
      tmpmin = movmin(data_in.trial{i}, winsize, 2);                        % get all minimum values
      if mod(winsize, 2)                                                    % remove useless results from the edges
        tmpmin = tmpmin(:, (winsize/2 + 1):(end-winsize/2));
      else
        tmpmin = tmpmin(:, (winsize/2 + 1):(end-winsize/2 + 1));
      end

      tmpmax = movmax(data_in.trial{i}, winsize, 2);                        % get all maximum values
      if mod(winsize, 2)                                                    % remove useless results from the edges
        tmpmax = tmpmax(:, (winsize/2 + 1):(end-winsize/2));
      else
        tmpmax = tmpmax(:, (winsize/2 + 1):(end-winsize/2 + 1));
      end

      tmpdiffmax = abs(tmpmax - tmpmedian);                                 % estimate the differences between the maximum values and the median
      tmpdiffmin = abs(tmpmin - tmpmedian);                                 % estimate the differences between the minimum values and the median
      tmp = cat(3, tmpdiffmax, tmpdiffmin);                                 % select always the maximum absolute difference
      tmp = max(tmp, [], 3);

      artfctmap{i} = tmp > cfgT.artfctdef.threshold.mad*tmpmad;             % find all violations
      [channum, begnum] = find(artfctmap{i});                               % estimate pairs of channel numbers and begin numbers for each violation
      artfctmap{i} = [artfctmap{i} false(length(channel), winsize - 1)];    % extend artfctmap to trial size
      endnum = begnum + winsize - 1;                                        % estimate end numbers for each violation
      for j=1:1:length(channum)
        artfctmap{i}(channum(j), begnum(j):endnum(j)) = true;               % extend the violations in the map to the window size
      end
      if ~isempty(begnum)
        begnum = unique(begnum);                                            % select all unique violations
        begnum = begnum + data_in.sampleinfo(i,1) - 1;                      % convert relative sample number into an absolute one
        begnum(:,2) = begnum(:,1) + winsize - 1;
        artifact = [artifact; begnum];                                      %#ok<AGROW>  add results to the artifacts matrix
      end
    end
  end

  autoart.artfctdef     = cfgT.artfctdef;                                   % generate output data structure
  autoart.showcallinfo  = cfgT.showcallinfo;
  autoart.artfctdef.threshold.artifact  = artifact;
  autoart.artfctdef.threshold.artfctmap = artfctmap;
  autoart.artfctdef.threshold.sliding   = 'yes';

end