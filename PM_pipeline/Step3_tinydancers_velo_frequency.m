%% Step 4 Tiny Dancers: Group analyses for frequencies
clear 
clc
% load all frequency data
cd('D:\TinyDancers_Movement\fieldtrip_velo\PM_PROC\Step2c_timefrequency\3m');
FileNames = {dir(['*_freq' '.mat']).name} ;
for k=1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_3m{k}     = freq_baseline;
    avg_control_3m{k}      = freq_control;
    avg_highvoice_3m{k}    = freq_highpitch;
    avg_lowpitch_3m{k}     = freq_lowpitch;

    cfg = [];
    cfg.parameter    = 'powspctrm';
    cfg.operation    = '(x1-x2)/(x1+x2)'; %we subtract the two power spectra and then divide them by their sum - this normalizes the difference by the common activity

    avg_mus_3m{k} = ft_math(cfg, freq_baseline, freq_control);
end

cd('D:\TinyDancers_Movement\fieldtrip_velo\PM_PROC\Step2c_timefrequency\6m');
FileNames = {dir(['*_freq' '.mat']).name} ;
for k=1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_6m{k}     = freq_baseline;
    avg_control_6m{k}      = freq_control;
    avg_highvoice_6m{k}    = freq_highpitch;
    avg_lowpitch_6m{k}     = freq_lowpitch;

    cfg = [];
    cfg.parameter    = 'powspctrm';
    cfg.operation    = '(x1-x2)/(x1+x2)'; %we subtract the two power spectra and then divide them by their sum - this normalizes the difference by the common activity

    avg_mus_6m{k} = ft_math(cfg, freq_baseline, freq_control);
end

cd('D:\TinyDancers_Movement\fieldtrip_velo\PM_PROC\Step2c_timefrequency\12m');
FileNames = {dir(['*_freq' '.mat']).name} ;
for k=1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_12m{k}     = freq_baseline;
    avg_control_12m{k}      = freq_control;
    avg_highvoice_12m{k}    = freq_highpitch;
    avg_lowpitch_12m{k}     = freq_lowpitch;

    cfg = [];
    cfg.parameter    = 'powspctrm';
    cfg.operation    = '(x1-x2)/(x1+x2)'; %we subtract the two power spectra and then divide them by their sum - this normalizes the difference by the common activity

    avg_mus_12m{k} = ft_math(cfg,freq_baseline, freq_control);
end

% calculate grand average for each condition
cfg = [];
cfg.channel = {'all'};
cfg.parameter = 'powspctrm';
grandavg_baseline_3m = ft_freqgrandaverage(cfg, avg_baseline_3m{:});
grandavg_control_3m = ft_freqgrandaverage(cfg, avg_control_3m{:});
grandavg_highvoice_3m = ft_freqgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowpitch_3m = ft_freqgrandaverage(cfg, avg_lowpitch_3m{:});
grandavg_mus_3m = ft_freqgrandaverage(cfg, avg_mus_3m{:});

grandavg_baseline_6m = ft_freqgrandaverage(cfg, avg_baseline_6m{:});
grandavg_control_6m = ft_freqgrandaverage(cfg, avg_control_6m{:});
grandavg_highvoice_6m = ft_freqgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowpitch_6m = ft_freqgrandaverage(cfg, avg_lowpitch_6m{:});
grandavg_mus_6m = ft_freqgrandaverage(cfg, avg_mus_6m{:});

grandavg_baseline_12m = ft_freqgrandaverage(cfg, avg_baseline_12m{:});
grandavg_control_12m = ft_freqgrandaverage(cfg, avg_control_12m{:});
grandavg_highvoice_12m = ft_freqgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowpitch_12m = ft_freqgrandaverage(cfg, avg_lowpitch_12m{:});
grandavg_mus_12m = ft_freqgrandaverage(cfg, avg_mus_12m{:});


%% 
cfg = [];
% cfg.baseline     = [-0.1 0];
% cfg.baselinetype = 'absolute';
% cfg.showlabels   = 'yes';
cfg.channel      = {'all'};
cfg.layout       = 'layout_sing.mat';
% cfg.xlim         = [0 5];
% cfg.ylim         = [0 5];
cfg.parameter    = 'powspctrm';
cfg.zlim       = [100 4000];
ft_singleplotTFR(cfg,grandavg_baseline_3m);%,grandavg_lowpitch_3m,grandavg_highvoice_3m);% grandavg_baseline,grandavg_control,grandavg_lowbass,grandavg_highvoice);
ft_singleplotTFR(cfg,grandavg_control_3m);%,grandavg_lowpitch_3m,grandavg_highvoice_3m);% grandavg_baseline,grandavg_control,grandavg_lowbass,grandavg_highvoice);
ft_singleplotTFR(cfg,grandavg_baseline_6m);
ft_singleplotTFR(cfg,grandavg_control_6m);
ft_singleplotTFR(cfg,grandavg_baseline_12m);
ft_singleplotTFR(cfg,grandavg_control_12m);


%%



cfg = [];
cfg.channel      = {'all'};
cfg.layout       = 'layout_sing.mat';
cfg.parameter    = 'powspctrm';

ft_singleplotTFR(cfg,grandavg_mus_3m);%,grandavg_lowpitch_3m,grandavg_highvoice_3m);% grandavg_baseline,grandavg_control,grandavg_lowbass,grandavg_highvoice);
ft_singleplotTFR(cfg,grandavg_mus_6m);%,grandavg_lowpitch_3m,grandavg_highvoice_3m);% grandavg_baseline,grandavg_control,grandavg_lowbass,grandavg_highvoice);
ft_singleplotTFR(cfg,grandavg_mus_12m);%,grandavg_lowpitch_3m,grandavg_highvoice_3m);% grandavg_baseline,grandavg_control,grandavg_lowbass,grandavg_highvoice);


%% stats
cfg                  = [];
% cfg.latency          = [0 0.444];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.05;      % alpha level of the permutation test
cfg.avgoverchan      = 'no';
cfg.tail             = 1;      % alpha level of the permutation test
cfg.clustertail      = 1;
cfg.numrandomization = 1000;
cfg.minnbchan        = 1;
neighbours(1).label = 'PM1';
neighbours(2).label = 'PM2';
neighbours(3).label = 'PM3';
neighbours(4).label = 'PM4';
neighbours(5).label = 'PM5';
neighbours(6).label = 'PM6';
neighbours(7).label = 'PM7';
neighbours(8).label = 'PM8';
neighbours(9).label = 'PM9';
neighbours(10).label = 'PM10';
% neighbours(1).neighblabel = {'PM2', 'PM5', 'PM6'}; % neighbours according to upper body, whole and lower body orga
% neighbours(2).neighblabel = {'PM1', 'PM5', 'PM6'};
% neighbours(3).neighblabel = {'PM8'};
% neighbours(4).neighblabel = {'PM7', 'PM9', 'PM10'};
% neighbours(5).neighblabel = {'PM2', 'PM1', 'PM6'};
% neighbours(6).neighblabel = {'PM2', 'PM1', 'PM5'};
% neighbours(7).neighblabel = {'PM4', 'PM9', 'PM10'};
% neighbours(8).neighblabel = {'PM3'};
% neighbours(9).neighblabel = {'PM7', 'PM4', 'PM10'};
% neighbours(10).neighblabel = {'PM7', 'PM4', 'PM9'};
neighbours(1).neighblabel = {'PM2','PM3','PM4', 'PM5', 'PM6','PM7', 'PM8', 'PM9', 'PM10'};% every PM orga can happen
neighbours(2).neighblabel = {'PM1','PM3','PM4', 'PM5', 'PM6','PM7', 'PM8', 'PM9', 'PM10'};
neighbours(3).neighblabel = {'PM1','PM2','PM4', 'PM5', 'PM6','PM7', 'PM8', 'PM9', 'PM10'};
neighbours(4).neighblabel = {'PM1','PM2','PM3', 'PM5', 'PM6','PM7', 'PM8', 'PM9', 'PM10'};
neighbours(5).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM6','PM7', 'PM8', 'PM9', 'PM10'};
neighbours(6).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM5','PM7', 'PM8', 'PM9', 'PM10'};
neighbours(7).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM5','PM6', 'PM8', 'PM9', 'PM10'};
neighbours(8).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM5','PM6', 'PM7', 'PM9', 'PM10'};
neighbours(9).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM5','PM6', 'PM8', 'PM8', 'PM10'};
neighbours(10).neighblabel = {'PM1','PM2','PM3', 'PM4', 'PM5','PM6', 'PM8', 'PM8', 'PM9'};
cfg.neighbours=neighbours;


Nsub=24;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;

stat_3m_mus = ft_freqstatistics(cfg, avg_baseline_3m{:}, avg_control_3m{:});
stat_3m_freq = ft_freqstatistics(cfg, avg_highvoice_3m{:}, avg_lowpitch_3m{:});


stat_12m_mus = ft_freqstatistics(cfg, avg_baseline_12m{:}, avg_control_12m{:});
stat_12m_freq = ft_freqstatistics(cfg, avg_lowpitch_12m{:}, avg_highvoice_12m{:});

Nsub=25;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;
stat_6m_mus = ft_freqstatistics(cfg, avg_baseline_6m{:}, avg_control_6m{:});
stat_6m_freq = ft_freqstatistics(cfg, avg_highvoice_6m{:}, avg_lowpitch_6m{:});


%%
freq=grandavg_baseline_ssep_3m.freq;
m3_ssep=[mean(grandavg_baseline_ssep_3m.powspctrm,1);mean(grandavg_control_ssep_3m.powspctrm,1)];
subplot(4,1,1) 
bar(freq,m3_ssep,'grouped');
ylim([-20000,50000]);

m6_ssep=[mean(grandavg_baseline_ssep_6m.powspctrm,1);mean(grandavg_control_ssep_6m.powspctrm,1)];
subplot(4,1,2) 
bar(freq,m6_ssep,'grouped');
ylim([-20000,75000]);

m12_ssep=[mean(grandavg_baseline_ssep_12m.powspctrm,1);mean(grandavg_control_ssep_12m.powspctrm,1)];
subplot(4,1,3) 
bar(freq,m12_ssep,'grouped');
ylim([-15000,40000]);

ma_ssep=[mean(grandavg_baseline_ssep.powspctrm,1);mean(grandavg_control_ssep.powspctrm,1)];
subplot(4,1,4) 
bar(freq,ma_ssep,'grouped');
ylim([-4000,8000]);

figure;
freq=grandavg_baseline_ssep_3m.freq;
m3_ssep=[mean(grandavg_lowbass_ssep_3m.powspctrm,1);mean(grandavg_highvoice_ssep_3m.powspctrm,1)];
subplot(4,1,1) 
bar(freq,m3_ssep,'grouped');
ylim([-20000,50000]);

m6_ssep=[mean(grandavg_lowbass_ssep_6m.powspctrm,1);mean(grandavg_highvoice_ssep_6m.powspctrm,1)];
subplot(4,1,2) 
bar(freq,m6_ssep,'grouped');
ylim([-30000,75000]);

m12_ssep=[mean(grandavg_lowbass_ssep_12m.powspctrm,1);mean(grandavg_highvoice_ssep_12m.powspctrm,1)];
subplot(4,1,3) 
bar(freq,m12_ssep,'grouped');
ylim([-15000,40000]);

ma_ssep=[mean(grandavg_lowbass_ssep.powspctrm,1);mean(grandavg_highvoice_ssep.powspctrm,1)];
subplot(4,1,4) 
bar(freq,ma_ssep,'grouped');
ylim([-4000,8000]);