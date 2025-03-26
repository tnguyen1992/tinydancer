function bad_channels = RFT_clean_flatlines(signal,max_flatline_duration,max_allowed_jitter)
% Detects (near-) flat-lined channels in your EEG data.
%
% ADAPTED BY NPA LAB to only return the list of bad channels, not removing
% them (as it was done in EEGlab and we work with ft structures)
%
% In:
%   Signal : continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%            with a 0.5Hz - 2.0Hz transition band)
%
%   MaxFlatlineDuration : Maximum tolerated flatline duration. In seconds. If a channel has a longer
%                         flatline than this, it will be considered abnormal. Default: 5
%
%   MaxAllowedJitter : Maximum tolerated jitter during flatlines. As a multiple of epsilon.
%                      Default: 20
%
% Out:
%   bad_channels : Cell of labels of the bad channels (then to be removed)
%
% Examples:
%   % use with defaults
%   eeg = clean_flatlines(eeg);
%
% Original authors               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-30
% Copyright (C) Christian Kothe, SCCN, 2012, ckothe@ucsd.edu 


if ~exist('max_flatline_duration','var') || isempty(max_flatline_duration) max_flatline_duration = 5; end
if ~exist('max_allowed_jitter','var') || isempty(max_allowed_jitter) max_allowed_jitter = 20; end

% Detect bad channels
bad_channels = {};
for c=1:signal.nbchan
    zero_intervals = reshape(find(diff([false abs(diff(signal.data(c,:)))<(max_allowed_jitter*eps) false])),2,[])';
    if max(zero_intervals(:,2) - zero_intervals(:,1)) > max_flatline_duration*signal.srate      % [original eeglab code]
        bad_channels{end+1} = signal.chanlocs(c).labels;     % [npa] finds label of the bad channel
    end
end
