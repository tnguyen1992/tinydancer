%% Step 4 Tiny Dancers: Group analyses for autocorrelation

%% Load data
clear
clc
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor\3m\');
respdata_mus=[];
respdata_freq=[];

FileNames = {dir(['*' '.mat']).name} ;
for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_3m{k}     = data_bas_autocor;
    avg_control_3m{k}      = data_con_autocor;
    clear data_bas_autocor data_con_autocor
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor_freq\3m\');
FileNames = {dir(['*' '.mat']).name} ;
for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_3m{k}    = data_hvo_autocor;
    avg_lowbass_3m{k}      = data_lba_autocor;

    clear data_hvo_autocor data_lba_autocor 
end


cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor\6m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})

    avg_baseline_6m{k}     = data_bas_autocor;
    avg_control_6m{k}      = data_con_autocor;
    
    clear data_bas_autocor data_con_autocor
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor_freq\6m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_6m{k}    = data_hvo_autocor;
    avg_lowbass_6m{k}      = data_lba_autocor;

    clear data_hvo_autocor data_lba_autocor
end




cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor\12m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})

    avg_baseline_12m{k}     = data_bas_autocor;
    avg_control_12m{k}      = data_con_autocor;

    clear data_bas_autocor data_con_autocor 
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor_freq\12m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_12m{k}    = data_hvo_autocor;
    avg_lowbass_12m{k}      = data_lba_autocor;

    clear data_hvo_autocor data_lba_autocor
end

%% grand average

cfg = [];
cfg.channel = {'all'};
% cfg.latency = [-0.1 1.0];
cfg.parameter = 'avg';
% cfg.nanmean        = 'yes';
grandavg_control_3m = ft_timelockgrandaverage_mod(cfg, avg_control_3m{:});
grandavg_control_6m = ft_timelockgrandaverage_mod(cfg, avg_control_6m{:});
grandavg_control_12m = ft_timelockgrandaverage_mod(cfg, avg_control_12m{:});

grandavg_baseline_12m = ft_timelockgrandaverage_mod(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m = ft_timelockgrandaverage_mod(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_timelockgrandaverage_mod(cfg, avg_lowbass_12m{:});

grandavg_highvoice_6m = ft_timelockgrandaverage_mod(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_timelockgrandaverage_mod(cfg, avg_lowbass_6m{:});
grandavg_baseline_6m = ft_timelockgrandaverage_mod(cfg, avg_baseline_6m{:});


grandavg_highvoice_3m = ft_timelockgrandaverage_mod(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_timelockgrandaverage_mod(cfg, avg_lowbass_3m{:});
grandavg_baseline_3m = ft_timelockgrandaverage_mod(cfg, avg_baseline_3m{:});

%% plot 
figure;
for k = 1:10
    subplot(2,5,k)
    time_x=grandavg_baseline_3m.time(1,:);
    y=grandavg_baseline_3m.avg(k,:);
    y_upper=y+(std(y)/sqrt(24));
    y_lower=y-(std(y)/sqrt(24));
    fill([time_x, fliplr(time_x)], [y_upper, fliplr(y_lower)],'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on
    plot(time_x,y)
    %     ylim([20 150]);
    title(sprintf("3m Movement to Tones in PM %d", k))
    y2=grandavg_control_3m.avg(k,:);
    y2_upper=y2+(std(y2)/sqrt(24));
    y2_lower=y2-(std(y2)/sqrt(24));
    fill([time_x, fliplr(time_x)], [y2_upper, fliplr(y2_lower)],'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_x,y2)
    %     ylim([20 150]);
    hold off
end

figure;
for k = 1:10
    subplot(2,5,k)
    time_x=grandavg_baseline_6m.time(1,:);
    y=grandavg_baseline_6m.avg(k,:);
    y_upper=y+(std(y)/sqrt(25));
    y_lower=y-(std(y)/sqrt(25));
    fill([time_x, fliplr(time_x)], [y_upper, fliplr(y_lower)],'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on
    plot(time_x,y)
    %     ylim([20 150]);
    title(sprintf("6m Movement to Tones in PM %d", k))
    y2=grandavg_control_6m.avg(k,:);
    y2_upper=y2+(std(y2)/sqrt(25));
    y2_lower=y2-(std(y2)/sqrt(25));
    fill([time_x, fliplr(time_x)], [y2_upper, fliplr(y2_lower)],'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_x,y2)
    %     ylim([20 150]);
    hold off
end
figure;
for k = 1:10
    subplot(2,5,k)
    time_x=grandavg_baseline_12m.time(1,:);
    y=grandavg_baseline_12m.avg(k,:);
    y_upper=y+(std(y)/sqrt(24));
    y_lower=y-(std(y)/sqrt(24));
    fill([time_x, fliplr(time_x)], [y_upper, fliplr(y_lower)],'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on
    plot(time_x,y)
    %     ylim([20 150]);
    title(sprintf("12m Movement to Tones in PM %d", k))
    y2=grandavg_control_12m.avg(k,:);
    y2_upper=y2+(std(y2)/sqrt(24));
    y2_lower=y2-(std(y2)/sqrt(24));
    fill([time_x, fliplr(time_x)], [y2_upper, fliplr(y2_lower)],'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_x,y2)
    %     ylim([20 150]);
    hold off
end

%% plot ERPstyle
cfg = [];
cfg.layout = 'layout_sing.mat';
cfg.channel = {'all'};

% cfg.xlim = [-0.1 0.6];
cfg.ylim = [0 0.1];
% ft_singleplotER(cfg,avg_baseline_3m{:})
% ft_singleplotER(cfg,avg_baseline_6m{:})
% ft_singleplotER(cfg,avg_baseline_12m{:})


[PKS,LOCS]=findpeaks(grandavg_baseline_3m.avg(1,:),25)
[PKS,LOCS]=findpeaks(grandavg_baseline_6m.avg(1,:),25)
[PKS,LOCS]=findpeaks(grandavg_baseline_12m.avg(1,:),25)

[PKS,LOCS]=findpeaks(grandavg_control_3m.avg(1,:),25)
[PKS,LOCS]=findpeaks(grandavg_control_6m.avg(1,:),25)
[PKS,LOCS]=findpeaks(grandavg_control_12m.avg(1,:),25)

ft_singleplotER(cfg,grandavg_baseline_3m,grandavg_control_3m)
ft_singleplotER(cfg,grandavg_baseline_6m,grandavg_control_6m)
ft_singleplotER(cfg,grandavg_baseline_12m,grandavg_control_12m)

%%
ft_singleplotER(cfg,grandavg_lowbass_3m, grandavg_highvoice_3m)
ft_singleplotER(cfg,grandavg_lowbass_6m, grandavg_highvoice_6m)
ft_singleplotER(cfg,grandavg_lowbass_12m, grandavg_highvoice_12m)

%%

cfg                  = [];
cfg.latency          = [0.400 0.520];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
cfg.correctm         = 'fdr';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.05;      % alpha level of the permutation test
% cfg.channel      = {'PM1', 'PM4'};
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.tail             = 1;      % alpha level of the permutation test
% cfg.clustertail      = 1;
cfg.numrandomization = 1000;
% cfg.minnbchan        = 1;
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

%24
Nsub=length(avg_highvoice_3m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;
stat_3m_freq = ft_timelockstatistics(cfg, avg_highvoice_3m{:}, avg_lowbass_3m{:});

Nsub=length(avg_lowbass_12m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;
stat_12m_freq = ft_timelockstatistics(cfg, avg_lowbass_12m{:}, avg_highvoice_12m{:});

Nsub=length(avg_baseline_12m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;

stat_12m_mus = ft_timelockstatistics(cfg, avg_baseline_12m{:}, avg_control_12m{:});

Nsub=length(avg_baseline_3m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;

stat_3m_mus = ft_timelockstatistics(cfg, avg_baseline_3m{:}, avg_control_3m{:});

%25

Nsub=length(avg_baseline_6m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;

stat_6m_mus = ft_timelockstatistics(cfg, avg_baseline_6m{:}, avg_control_6m{:});

Nsub=length(avg_highvoice_6m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;
stat_6m_freq = ft_timelockstatistics(cfg, avg_highvoice_6m{:}, avg_lowbass_6m{:});







