%% Step 1 Tiny Dancers: Load EEG and preprocess
% Steps are Bandpass filtering, Bad channel rejection, ASR,
% Automatic ICA, Bad channel rejection part 2 and interpolation, followed
% by rereferencing to the average
% Trinh Nov 22

%% initialise fieldtrip
clc
clear all
close all

addpath('E:/MUSICOM_TinyDancer/Adults/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
addpath(genpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/Giac_ToolBox/'))
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/clean_rawdata2.7')
addpath 'D:/MUSICOM_EEG/onsets'
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/functions')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/external/eeglab/')
ft_defaults;
%% Load, Epoching and BP filter
% go to data folder
cd(['E:/MUSICOM_TinyDancer/Adults/EEG/']);
% load names of all files
FileNames      = {dir(['*' '.bdf']).name} ;
group='adults';
%% start loop
for k = 27:1:length(FileNames)
    % get baby ID
    baby_name      = erase(FileNames{k},'.bdf');

    % load silence condition
    data_silence   = Trinh_loadConditionData(FileNames{k},'STATUS',{1},3,13);
    % load auditory conditions
    data_auditory  = Trinh_loadConditionData(FileNames{k},'STATUS',{2,3,4,5,6,7,8,9},3,24);

    % append silence and auditory conditions
    cfg            = [];
    data_all       = ft_appenddata(cfg, data_silence, data_auditory);

    % bandpass filter data
    data_hpfilt      = Trinh_filterEEG(data_all,'hp',[0.3]);
    data_filt        = Trinh_filterEEG(data_hpfilt,'lp',[30]);


    %% extract VEOG for ICA and rename channels
    cfg            = [];
    cfg.channel    = {'EXG1'};
    data_veog      = ft_selectdata(cfg,data_all); % saves the veog for later
    %     data_woveog    = Giac_removeChannels(data_all, {'EXG1'} );
    data_woveog    = Giac_removeChannels(data_filt, {'EXG1','EXG4', ...
        'EXG5','EXG6','EXG7','EXG8', 'Status'} );

    % save channel configuration for interpolation
    [channel]      = ft_channelselection({'all','-EXG3', '-EXG3'},data_woveog, 'EEG');

    data_woveog.label{65, 1}='TP9';
    data_woveog.label{66, 1}='TP10';




    %% Exclude bad channels and trials
    % check for flatlines and noisy channels according to SD, mean and peak to
    % peak
    steps_to_do    = {'flat','corr-noise'};
    criteria       = {[5],[0.1,20]};
    flat_channels  = RFT_detect_FlatCorrNoise_electrodes(data_woveog,'Trinh_tinydancer_eeglab_template.mat',steps_to_do,criteria);
    data_flat      = Giac_removeChannels(data_woveog, flat_channels.flat );
    try
        data_flat      = Giac_removeChannels(data_flat, flat_channels.corrNoise );
    end

    % get noisy bad channels according to SD,mean and peak to peak info
    bad_channels   = Giac_EEG_CatchNoisyElectrodes(data_flat, 'all', 3, 'recursive' );
    data_badchan   = Giac_removeChannels(data_flat, bad_channels );


    %% reref
    try
        data_reref = Trinh_reref(data_badchan,[] , {'TP9', 'TP10'});%{'all','-VEOG'});
    catch
        fprintf('TP10 missing!')
    end

    data_reref = Giac_removeChannels(data_reref, {'TP9'} );

    % visualise data
     Trinh_Visualise(data_reref,'layout_sing.mat',[-30 30])
    %% Artefact correction
    % ASR seems better thus far,
    % wavthresh introduces some weird slow waves and seems to extract a lot of
    % neural data or saturates channels unneccesarily
    data_asr       = RFT_clean_asr_combined_trials_ftstruct(data_reref,'Trinh_tinydancer_eeglab_template.mat', 'layout_sing.mat', 5,[],[],[]);

    % Quick visualization
%     Trinh_Visualise(data_asr,'layout_sing.mat',[-30 30])

    %% ICA
    % cfg          = [];
    % data_prepICA = ft_appenddata(cfg, data_asr, data_veog);
    [data_ICA,~] = RFT_IClabel(data_asr,'Trinh_tinydancer_eeglab_template.mat',30,[0 0;0 0; 0.5 1; 0 0; 0 0; 0 0; 0 0]);

    Trinh_Visualise(data_ICA,'layout_sing.mat',[-30 30])
%     close all
    %% check bad channels again and interpolate if still bad

    bad_channels      = Giac_EEG_CatchNoisyElectrodes(data_ICA, 'all', 3, 'recursive' );
    data_prepint      = Giac_removeChannels( data_ICA, bad_channels );

    % get missing channels by subtracting new channels from old channels
%     addpath 'D:\MUSICOM_EEG'
    load('Trinh_tinydancer_channels.mat');
    chanIDs                     = [];
    chanIDs                     = setdiff(channel,data_prepint.label ) ;
    badchannel_overview{k,1}    = chanIDs;

    % interpolate bad channels
    data_int                    = Trinh_interpolateChannels(data_prepint,chanIDs,'layout_sing_neighb.mat','spline');

    %% save data
    save(strcat('D:\MUSICOM_EEG\PROC\Step1_processeddatamastoidreref/',group,'/',baby_name,'_PROC'),'data_int');
end





