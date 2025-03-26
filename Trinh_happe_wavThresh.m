function ft_struct = Trinh_happe_wavThresh(ft_struct,wavLvl,wavFam,ThresholdRule)
% function uses the wavelet thresholding algorithm from HAPPILEE
% needed input: ft_struct from ft_preprocessing
% wavelevel 10 for sampling rate higher than 500, 9 for 250-500, 8 for
% lower than 250
% Set the wavelet family depending on user input wavFam = 'bior5.5'(default) 
% otherwise: wavFam = 'coif4' ;
% Set the threshold rule depending on user input (when applicable) and
% paradigm: 'Soft' or 'Hard'
% 12/10/22 Trinh Nguyen (IIT)

% get number of channels and reshape ft to EEGLAB format
[nChan] = size(ft_struct.trial{1,1},1);
all_trials = cat(2,ft_struct.trial{:});

% Use wavelet thresholding to determine artifacts in the data.
artifacts = wdenoise(reshape(all_trials, size(all_trials, 1), [])', wavLvl, ...
    'Wavelet', wavFam, 'DenoisingMethod', 'Bayes', 'ThresholdRule', ...
    ThresholdRule, 'NoiseEstimate', 'LevelDependent')' ;
    
% substract artifacts from EEG
all_trials_post = all_trials - artifacts ;

% reformat EEG data structure back to ft (uneven trials allowed)
ft_struct.trial = mat2cell(all_trials_post,nChan,cell2mat(cellfun(@(x) size(x,2),ft_struct.trial,'UniformOutput',0)));
end