%% Frequency analyses for Tiny dancers study
clear all
clc

%% access the continous data
group={'3m','6m','12m','adults'};
for g=1:length(group)
    cd(['D:\MUSICOM_EEG\PROC\Step2_continuousdata\' group{g}])
    FileNames = {dir(['*_CONTIN' '.mat']).name} ;
    for k=1:1:length(FileNames)

        %% get baby ID and load data
        baby_name = erase(FileNames{k},'_CONTIN.mat');
        load(FileNames{k}) ;

        cfg                          = [];
        data_baseline_contin         = ft_timelockanalysis(cfg, data_baseline_contin);
        data_control_contin          = ft_timelockanalysis(cfg, data_control_contin);
        data_lowbass_contin          = ft_timelockanalysis(cfg, data_lowbass_contin);
        data_highvoice_contin        = ft_timelockanalysis(cfg, data_highvoice_contin);

        %% run frequency analysis
        cfg             = [];
%         cfg.output      = 'pow';
        cfg.output      = 'fourier';
        cfg.channel     = {'all'};
        cfg.method      = 'mtmfft';
        cfg.taper       = 'hilbert';
%                 cfg.foilim      = [1 5];
        cfg.foi      = 1:0.125:10;
%         cfg.foi      = 1:0.5:15;

        %         cfg.polyremoval = 1;
        %         cfg.tapsmofrq   = 4;
        %         cfg.pad         ='nextpow2';

        freq_baseline   = ft_freqanalysis(cfg, data_baseline_contin);
        %         freq_baseline   = ft_freqanalysis(cfg, data_out);

        %         plot(freq_baseline.freq, squeeze(abs(freq_baseline.fourierspctrm(1,1,:))));
        % plot(freq_baseline.freq, freq_baseline.powspctrm(1,:))
        %                 hold on
        %         %
        %         cfg.taper       = 'boxcar';
        %         freq_baseline_h   = ft_freqanalysis(cfg, data_baseline_contin);
        % %         plot(freq_baseline_h.freq, squeeze(abs(freq_baseline_h.fourierspctrm(1,1,:))));
        %         plot(freq_baseline_h.freq, freq_baseline_h.powspctrm(1,:))
        freq_control    = ft_freqanalysis(cfg, data_control_contin);

        freq_lowbass    = ft_freqanalysis(cfg, data_lowbass_contin);

        freq_highvoice  = ft_freqanalysis(cfg, data_highvoice_contin);

        cfg             = [];
        [freq_baseline] = ft_freqdescriptives(cfg, freq_baseline);
        [freq_control] = ft_freqdescriptives(cfg, freq_control);
        [freq_lowbass] = ft_freqdescriptives(cfg, freq_lowbass);
        [freq_highvoice] = ft_freqdescriptives(cfg, freq_highvoice);

        save(strcat('D:\MUSICOM_EEG\PROC\Step3_frequencydata\',group{g},'\',baby_name,'_FREQ'),...
            'freq_baseline','freq_control','freq_lowbass','freq_highvoice');

    end
end
