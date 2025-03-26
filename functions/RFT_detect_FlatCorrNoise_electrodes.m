function [bad_channels] = RFT_detect_FlatCorrNoise_electrodes(ft_struct,eeglab_template_file,steps_to_do,criteria)
%%% From a fieldtrip structure (with N trials), concatenates all trials,
%%% create an eeglab structure with this 'single concatenated trial',
%%% Detects channels that are either 'flat-lined', 'line-noise +
%%% low-correlation', or both in your EEG data.
%%%
%%% Inputs:
%%% ft_struct                 Fieldtrip struct with the original data
%%% eeglab_template_file      Name of the file of your eeglab template
%%% step_to_do                Cell of strings stating the cleaning steps
%%% criteria                  Cell of criteria for each cleaning step
%%%
%%% Output:
%%% bad_channels              Struct (fields 'flat','corr-noise') with cells of labels of the bad channels
%%% NPA - Oct 2022 %%%

assert( (exist('steps_to_do','var') && ~isempty(steps_to_do)) , 'cleaning not possible if no steps specified');
assert( (exist('criteria','var') && ~isempty(criteria)) , 'cleaning not possible if no criteria specified');


%%%% CONCATENATE ALL YOUR TRIALS %%%%
all_trials = cat(2,ft_struct.trial{:});


%%%%    CONVERT FIELDTRIP STRUCT TO EEGLAB STRUCT     %%%%
%%%% (! might need to be validated and/or improved !) %%%%
%%%% (inspired by EEGlab func. 'fieldtrip2eeglab.m')  %%%%

%%% load template to create eeglab structure %%%
eeglab_template = load(eeglab_template_file);
try
    eeglab_struct = eeglab_template.eeglab_template;
catch
    eeglab_struct = eeglab_template.tinydancer_eeglab_template;
end

%%% update the data field (Nchan x Ntime) %%%
eeglab_struct.data = all_trials;

%%% update other fields of interest %%%
[eeglab_struct.nbchan , eeglab_struct.pnts] = size(all_trials);         % number of channels , number of time points
eeglab_struct.trials = 1;           % because 1 big trial (of concatenated data)
eeglab_struct.srate = ft_struct.fsample;         % important if you have a sampling rate ≠ 2048Hz

% The following part is maybe not necessary (it updates the time
% information, saying "tStart = tStart of the original data in trial 1" and
% then deducting tStop based on srate and the num of time points concatenated
% (eg, you concatenated 3 trials [-1s, 10s], it's gonna create [-1s,32s]))
eeglab_struct.xmin = ft_struct.time{1}(1);
eeglab_struct.xmax =  eeglab_struct.xmin + (eeglab_struct.pnts-1)/eeglab_struct.srate;
eeglab_struct.times = eeglab_struct.xmin : 1/eeglab_struct.srate : eeglab_struct.xmax;

%%% update chanlocs if you have removed channels %%%
chans_removed_in_ft=[];
for ch=1:length(eeglab_struct.chanlocs)
    if ~ismember(eeglab_struct.chanlocs(ch).labels,ft_struct.label)   chans_removed_in_ft(end+1)=ch;        end
end
eeglab_struct.chanlocs(chans_removed_in_ft) = [];


%%%% CALL THE CLEANING FUNCTIONS %%%%

%%% Flat line detection
if any(strcmp('flat',steps_to_do))
    disp('Detecting flat line...')
    idxFlat = find(strcmp('flat',steps_to_do));
    flatline_crit = criteria{idxFlat};
    bad_channels.flat = RFT_clean_flatlines(eeglab_struct,flatline_crit);

    if ~isempty(bad_channels.flat)
        disp('Bad channels detected (flat line): ');  disp(bad_channels.flat);
    else  disp('No bad channels detected (flat line)');
    end
end

%%% Noise and low-correlation detection
if any(strcmp('corr-noise',steps_to_do))
    disp('Detecting line-noise and low-correlation with neighbours...')
    idxCorrNoise = find(strcmp('corr-noise',steps_to_do));
    corr_crit = criteria{idxCorrNoise}(1);
    noise_crit = criteria{idxCorrNoise}(2);
    [~,noisy_channels] = clean_channels(eeglab_struct,corr_crit,noise_crit);

    noisy_struct=eeglab_struct.chanlocs(noisy_channels);
    noisy_cell=struct2cell(noisy_struct);
    if ~isempty(noisy_cell)
        bad_channels.corrNoise={noisy_cell{1,1,:}};
        disp('Bad channels detected (corr-noise): ');  disp(bad_channels.corrNoise);
    else  disp('No bad channels detected (corr-noise)');
        bad_channels.corrNoise=[];
    end


end

end
