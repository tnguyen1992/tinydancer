%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear all
clc

addpath 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\functions'
addpath(genpath( 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\toolboxes/Giac_ToolBox/'))

addpath 'D:\MUSICOM_EEG/onsets'
addpath('D:\MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
ft_defaults;
%% access preprocessed files
group={'3m','6m','12m'};
for g=1:length(group)
    cd(['D:\TinyDancers_Movement\fieldtrip_velo\POS\' group{g}]);
    FileNames = {dir(['*' '.mat']).name} ;

    for k=1:1:length(FileNames)

        % load data
        load(FileNames{k});
        baby_name       = erase(FileNames{k},'.mat');

        % rectify the velocity
        cfg=[];
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.5;
        data_baseline = ft_preprocessing(cfg,data_baseline);
        data_control = ft_preprocessing(cfg,data_control);
        data_highvoice = ft_preprocessing(cfg,data_highvoice);
        data_lowbass = ft_preprocessing(cfg,data_lowbass);

%         % smoothing data
        data_baseline = Giac_MovingMean( data_baseline, 3 );
        data_control = Giac_MovingMean( data_control, 3 );
        data_highvoice = Giac_MovingMean( data_highvoice, 3 );
        data_lowbass = Giac_MovingMean( data_lowbass, 3 );

        %% window
        cfg = [];
        cfg.length  = 7;
        cfg.overlap = 0;

        data_baseline_win   = ft_redefinetrial(cfg, data_baseline);
        data_control_win    = ft_redefinetrial(cfg, data_control);
        data_lowbass_win    = ft_redefinetrial(cfg, data_lowbass);
        data_highvoice_win  = ft_redefinetrial(cfg, data_highvoice);

        %% FFT
        cfg             = [];
        cfg.output      = 'fourier';
        cfg.channel     = {'all'};
        cfg.method      = 'mtmfft';
        cfg.taper       = 'hilbert';
        cfg.foi         = [1.500:0.125:5];

        freq_baseline   = ft_freqanalysis(cfg, data_baseline_win);
        freq_control    = ft_freqanalysis(cfg, data_control_win);
        freq_highpitch  = ft_freqanalysis(cfg, data_highvoice_win);
        freq_lowpitch   = ft_freqanalysis(cfg, data_lowbass_win);

        [freq_baseline] = ft_freqdescriptives(cfg, freq_baseline);
        [freq_control]  = ft_freqdescriptives(cfg, freq_control);
        [freq_lowpitch] = ft_freqdescriptives(cfg, freq_lowpitch);
        [freq_highpitch]= ft_freqdescriptives(cfg, freq_highpitch);

        save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\POS_PROC\Step2b_frequency\',group{g},'\',baby_name,'_freq'),...
            'freq_baseline','freq_control','freq_lowpitch','freq_highpitch');
    end
end

