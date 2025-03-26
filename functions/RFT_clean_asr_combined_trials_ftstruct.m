function [ft_struct] = RFT_clean_asr_combined_trials_ftstruct(ft_struct,eeglab_template_file, ft_template_file, cutoff,windowlen,ref_maxbadchannels,ref_tolerances)
%%% From a fieldtrip structure (with N trials), concatenates all trials,   
%%% create an eeglab structure with this 'single concatenated trial',       
%%% cleans the data with ASR, based on cutoff and windowlen                 
%%%                                                                         
%%% Inputs:
%%% ft_struct                 Fieldtrip struct with the original data
%%% eeglab_template_file      Name of the file of your eeglab template
%%% cutoff                    Cutoff parameter for ASR datacleaning
%%% windowlen                 Window length parameter for ASR datacleaning
%%% ref_tolerances            SD of outliers in power determined for the
%%%                           reference data
%%% 
%%% Output:
%%% ft_struct                 Fieldtrip struct with the asr-cleaned data
%%% Atesh & Fefe & Trinh & Robs - Oct 2022 %%%

if nargin <4
    cutoff = [];
end

if nargin <5
    windowlen = [];
end

if nargin <6
    ref_maxbadchannels = [];
end

if nargin <7
    ref_tolerances = [];
end

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



%%%% CALL ASR %%%%
usegpu = true;        % decide whether you run ASR on GPU or not
%%% RB: modified clean_asr function to get the cleaned reference data as
%%% output
[asr_cleaned_eeglab_struct, ref_section] = RFT_clean_asr(eeglab_struct,cutoff,windowlen,[],[],ref_maxbadchannels,ref_tolerances,[],usegpu);  % asr
asr_cleaned_combined_trials = asr_cleaned_eeglab_struct.data;       % take the cleaned concatenated data from eeglab struct


%%%% SEGMENT THE DATA BACK IN TRIALS %%%%
[nChan,nTimePoints] = size(ft_struct.trial{1,1});
ft_struct.trial = mat2cell(asr_cleaned_combined_trials,nChan,cell2mat(cellfun(@(x) size(x,2),ft_struct.trial,'UniformOutput',0)));      

% %%%% RB: visualise clean reference data for ASR with whole data 
% trial = ref_section.data; 
% ft_dummy = struct;
% ft_dummy.trial{1} = trial;
% ft_dummy.fsample = ft_struct.fsample;
% ft_dummy.time{1} = linspace(0,size(trial,2)/ft_dummy.fsample, size(trial,2));
% ft_dummy. label = ft_struct.label;
% [~]= RFT_EEG_Visualise( ft_dummy, {'all'}, [-50 50], ft_template_file);
% trial = all_trials; 
% ft_dummy = struct;
% ft_dummy.trial{1} = trial;
% ft_dummy.fsample = ft_struct.fsample;
% ft_dummy.time{1} = linspace(0,size(trial,2)/ft_dummy.fsample, size(trial,2));
% ft_dummy. label = ft_struct.label;
% [~] = RFT_EEG_Visualise( ft_dummy, {'all'}, [-50 50], ft_template_file);
% trial = asr_cleaned_combined_trials; 
% ft_dummy = struct;
% ft_dummy.trial{1} = trial;
% ft_dummy.fsample = ft_struct.fsample;
% ft_dummy.time{1} = linspace(0,size(trial,2)/ft_dummy.fsample, size(trial,2));
% ft_dummy. label = ft_struct.label;
% [~] = RFT_EEG_Visualise( ft_dummy, {'all'}, [-50 50], ft_template_file);

end
