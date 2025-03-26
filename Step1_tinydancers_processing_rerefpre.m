%% Step 1 Tiny Dancers: Load EEG and preprocess
% Steps are Bandpass filtering, Bad channel rejection, ASR,
% Automatic ICA, Bad channel rejection part 2 and interpolation, followed
% by rereferencing to the average
% Trinh Nov 22

%% initialise fieldtrip
clc
clear all
close all

addpath('D:/MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
addpath(genpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/Giac_ToolBox/'))
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/clean_rawdata2.7')
addpath 'D:/MUSICOM_EEG/onsets'
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/functions')
% addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/external/eeglab/')
ft_defaults;
%% Load, Epoching and BP filter
group={'3m','6m','12m'};
countcomps=[];
for g=1:length(group)
    % go to data folder
    cd(['D:\MUSICOM_EEG\EEG/' group{g}]);
    % load names of all files
    FileNames      = {dir(['*' '.vhdr']).name} ;

    %% start loop
    for k = 1:1:length(FileNames)
        % get baby ID
        baby_name      = erase(FileNames{k},'.vhdr');

        % load silence condition
        data_silence   = Trinh_loadConditionData(FileNames{k},'Stimulus',{'S  1'},3,13);
        % load auditory conditions
        data_auditory  = Trinh_loadConditionData(FileNames{k},'Stimulus',{'S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9'},3,24);

        % append silence and auditory conditions
        cfg            = [];
        data_all       = ft_appenddata(cfg, data_silence, data_auditory);

        % bandpass filter data
        data_filt      = Trinh_filterEEG(data_all,'bp',[0.3 30]);

        % extract VEOG for ICA
        cfg            = [];
        cfg.channel    = {'VEOG'};
        data_veog      = ft_selectdata(cfg,data_all); % saves the veog for later
        data_woveog    = Giac_removeChannels(data_filt, {'VEOG'} );

        % save channel configuration for interpolation
        % [channel]      = ft_channelselection({'all','-VEOG'},data_all, 'EEG');


        %% Exclude bad channels and trials
        % check for flatlines and noisy channels according to SD, mean and peak to
        % peak
        steps_to_do    = {'flat','corr-noise'};
        criteria       = {[5],[0.1,20]};
        flat_channels  = RFT_detect_FlatCorrNoise_electrodes(data_woveog,'Trinh_MUSICOM_eeglab_template.mat',steps_to_do,criteria);
        data_flat      = Giac_removeChannels(data_woveog, flat_channels.flat );
        data_flat      = Giac_removeChannels(data_flat, flat_channels.corrNoise );

        % get noisy bad channels according to SD,mean and peak to peak info
        bad_channels   = Giac_EEG_CatchNoisyElectrodes(data_flat, 'all', 2.75, 'recursive' );
        data_badchan   = Giac_removeChannels(data_flat, bad_channels );


        %% reref
        try
            data_reref = Trinh_reref(data_badchan, 'TP9', {'TP9', 'TP10'});%{'all','-VEOG'});
            data_reref = Giac_removeChannels(data_reref, {'TP9'} );

            % visualise data
            %         Trinh_Visualise(data_reref,'layout_sing.mat',[-30 30])
            %% Artefact correction
            % ASR seems better thus far,
            % wavthresh introduces some weird slow waves and seems to extract a lot of
            % neural data or saturates channels unneccesarily
            data_asr       = RFT_clean_asr_combined_trials_ftstruct(data_reref,'Trinh_MUSICOM_eeglab_template.mat', 'layout_sing.mat', 5,[],[],[]);

            % Quick visualization
            %         Trinh_Visualise(data_asr,'layout_sing.mat',[-30 30])

            %% ICA
            % cfg          = [];
            % data_prepICA = ft_appenddata(cfg, data_asr, data_veog);
            [data_ICA,rejected_comps] = RFT_IClabel(data_asr,'Trinh_MUSICOM_eeglab_template.mat',size(data_asr.label,1)-1,[0 0;0 0; 0.5 1; 0 0; 0 0; 0 0; 0 0]);
            countcomps = [countcomps;length(rejected_comps)];
            % Trinh_Visualise(data_ICA,'layout_sing.mat',[-30 30])
            close all
            %% check bad channels again and interpolate if still bad

            bad_channels      = Giac_EEG_CatchNoisyElectrodes(data_ICA, 'all', 2.75, 'recursive' );
            data_prepint      = Giac_removeChannels( data_ICA, bad_channels );

            % get missing channels by subtracting new channels from old channels
            addpath 'D:\MUSICOM_EEG'
            load('MUSICOM_channels.mat');
            chanIDs                     = [];
            chanIDs                     = setdiff(channel,data_prepint.label ) ;
            badchannel_overview{k,1}    = chanIDs;

            % interpolate bad channels
            data_int                    = Trinh_interpolateChannels(data_prepint,chanIDs,'layout_sing_neighb.mat','spline');

            %% save data
            save(strcat('D:\MUSICOM_EEG\PROC\Step1_processeddatamastoidreref/',group{g},'/',baby_name,'_PROC'),'data_int');
        end
    end
end





