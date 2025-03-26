%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear all
clc

addpath 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\functions'
addpath 'D:\MUSICOM_EEG/onsets'
addpath('D:\MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
ft_defaults;
%% access preprocessed files
group={'3m','6m','12m'};
for g=1:length(group)
    cd(['D:\TinyDancers_Movement\fieldtrip_velo\PM\' group{g}]);
    FileNames = {dir(['*' '.mat']).name} ;

    for k=1:1:length(FileNames)

        % load data
        load(FileNames{k});
        baby_name       = erase(FileNames{k},'.mat');

        % rectify the velocity
%         cfg=[];
%         cfg.rectify = 'yes';
%         data_baseline = ft_preprocessing(cfg,data_baseline);
%         data_control = ft_preprocessing(cfg,data_control);
%         data_highvoice = ft_preprocessing(cfg,data_highvoice);
%         data_lowbass = ft_preprocessing(cfg,data_lowbass);

        % smoothing data
        data_baseline = Giac_MovingMean( data_baseline, 3 );
        data_control = Giac_MovingMean( data_control, 3 );
        data_highvoice = Giac_MovingMean( data_highvoice, 3 );
        data_lowbass = Giac_MovingMean( data_lowbass, 3 );

        data_baseline = Giac_MovingMean( data_baseline, 3 );
        data_control = Giac_MovingMean( data_control, 3 );
        data_highvoice = Giac_MovingMean( data_highvoice, 3 );
        data_lowbass = Giac_MovingMean( data_lowbass, 3 );
       

        %% window
%         cfg = [];
%         cfg.length  = 3.5;
%         cfg.overlap = 0.8;
% 
%         data_baseline_abs_superPM    = ft_redefinetrial(cfg, data_baseline);
%         data_control_abs_superPM    = ft_redefinetrial(cfg, data_control);
%         data_lowbass_abs_superPM    = ft_redefinetrial(cfg, data_lowbass);
%         data_highvoice_abs_superPM    = ft_redefinetrial(cfg, data_highvoice);

        %% FFT
        cfg             = [];
        cfg.output      = 'pow';
        cfg.channel     = {'all'};
        cfg.method      = 'mtmconvol';
        cfg.taper       = 'hanning';
        cfg.foi   = 0.5:0.05:5;
%         cfg.toi   = '80%';
        cfg.toi     = 0:0.04:21;
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*2;   % length of time window = 0.5 sec

        freq_baseline   = ft_freqanalysis(cfg, data_baseline);
        freq_control   = ft_freqanalysis(cfg, data_control);
        freq_highpitch   = ft_freqanalysis(cfg, data_highvoice);
        freq_lowpitch   = ft_freqanalysis(cfg, data_lowbass);

        [freq_baseline] = ft_freqdescriptives(cfg, freq_baseline);
        [freq_control] = ft_freqdescriptives(cfg, freq_control);
        [freq_lowpitch] = ft_freqdescriptives(cfg, freq_lowpitch);
        [freq_highpitch] = ft_freqdescriptives(cfg, freq_highpitch);

        save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\PM_PROC\Step2c_timefrequency\',group{g},'\',baby_name,'_freq'),...
            'freq_baseline','freq_control','freq_lowpitch','freq_highpitch');
    end
end

