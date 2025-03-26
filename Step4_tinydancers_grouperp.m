%% Step 4 Tiny Dancers: Group analyses for ERP
% load all ERP data
clear all
clc

%% Load data
cd('D:\MUSICOM_EEG\PROC\Step3_ERPs\3m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_3m{k}     = ERP_baseline;
    avg_control_3m{k}      = ERP_control;
    avg_highvoice_3m{k}    = ERP_highvoice;
    avg_lowbass_3m{k}      = ERP_lowbass;
%     avg_music_3m{k}        = ERP_music;
    clear ERP_baseline ERP_control ERP_highvoice ERP_lowbass
end

cd('D:\MUSICOM_EEG\PROC\Step3_ERPs\6m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_6m{k}     = ERP_baseline;
    avg_control_6m{k}      = ERP_control;
    avg_highvoice_6m{k}    = ERP_highvoice;
    avg_lowbass_6m{k}      = ERP_lowbass;
    avg_zero_6m{k}         = ERP_control;
    avg_zero_6m{k}.avg     = zeros(32,2000);
    clear ERP_baseline ERP_control ERP_highvoice ERP_lowbass
end

cd('D:\MUSICOM_EEG\PROC\Step3_ERPs\12m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_12m{k}     = ERP_baseline;
    avg_control_12m{k}      = ERP_control;
    avg_highvoice_12m{k}    = ERP_highvoice;
    avg_lowbass_12m{k}      = ERP_lowbass;
    clear ERP_baseline ERP_control ERP_highvoice ERP_lowbass
end

cd('D:\MUSICOM_EEG\PROC\Step3_ERPs\adults\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_ad{k}     = ERP_baseline;
    avg_control_ad{k}      = ERP_control;
    avg_highvoice_ad{k}    = ERP_highvoice;
    avg_lowbass_ad{k}      = ERP_lowbass;
    clear ERP_baseline ERP_control ERP_highvoice ERP_lowbass
end

avg_music = [avg_baseline_3m avg_baseline_6m avg_baseline_12m];
avg_control = [avg_control_3m avg_control_6m avg_control_12m];

%% calculate grand average for each condition
cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.1 0.5];
cfg.parameter = 'avg';

grandavg_control_3m = ft_timelockgrandaverage(cfg, avg_control_3m{:});
grandavg_control_6m = ft_timelockgrandaverage(cfg, avg_control_6m{:});
grandavg_control_12m = ft_timelockgrandaverage(cfg, avg_control_12m{:});
grandavg_control_ad = ft_timelockgrandaverage(cfg, avg_control_ad{:});

grandavg_baseline_12m = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});
grandavg_highvoice_6m = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_baseline_6m = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});
grandavg_highvoice_3m = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_baseline_3m = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});

grandavg_highvoice_ad = ft_timelockgrandaverage(cfg, avg_highvoice_ad{:});
grandavg_lowbass_ad = ft_timelockgrandaverage(cfg, avg_lowbass_ad{:});
grandavg_baseline_ad = ft_timelockgrandaverage(cfg, avg_baseline_ad{:});


%% plot ERPs
cfg = [];
cfg.layout = 'layout_sing.mat';
cfg.channel = {'C3','Cz','C4','FC3','FCz','FC4','F3','Fz','F4'};
cfg.xlim = [-0.1 0.6];
cfg.zlim = [-0.5 4];

ft_singleplotER(cfg,grandavg_baseline_3m,grandavg_control_3m)
% ylim([-2 4])
ft_singleplotER(cfg,grandavg_baseline_6m,grandavg_control_6m)
% ylim([-1.5 2.5])
ft_singleplotER(cfg,grandavg_baseline_12m,grandavg_control_12m)
% ylim([-1.5 2.5])
cfg.channel = {'C3','Cz','C4','FC3','FCz','FC4','F3','Fz','F4','CP3','CP4','P3',...
    'P4','Pz'};
ft_singleplotER(cfg,grandavg_baseline_ad,grandavg_control_ad)
ylim([-1 1.25])

% ft_singleplotER(cfg,grandavg_lowbass_3m, grandavg_highvoice_3m)
% ylim([-1.7 2.8])
% ft_singleplotER(cfg,grandavg_lowbass_6m, grandavg_highvoice_6m)
% ylim([-1.7 2.8])
% ft_singleplotER(cfg,grandavg_lowbass_12m, grandavg_highvoice_12m)
% ylim([-1.7 2.8])
% ft_singleplotER(cfg,grandavg_lowbass_ad, grandavg_highvoice_ad)
% ylim([-1 1.25])

%% plot Figure for paper
% Plot 3m data
figure;
time_x = grandavg_baseline_3m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_baseline_3m.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_baseline_3m.var([1:13 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_baseline_3m.var([1:13 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("3m ERP (Music vs Shuffled)");

% Average over all PMs for control condition
y_control = mean(grandavg_control_3m.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_control_3m.var([1:13 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_control_3m.var([1:13 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot 6m data
figure;
time_x = grandavg_baseline_6m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_baseline_6m.avg([1:18 24:27 ],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_baseline_6m.var([1:18 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_baseline_6m.var([1:18 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("6m ERP (Music vs Shuffled)");

% Average over all PMs for control condition
y_control = mean(grandavg_control_6m.avg([1:18 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_control_6m.var([1:18 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_control_6m.var([1:18 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot 12m data
figure;
time_x = grandavg_baseline_12m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_baseline_12m.avg([1:12 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_baseline_12m.var([1:12 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_baseline_12m.var([1:12 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("12m ERP (Music vs Shuffled)");

% Average over all PMs for control condition
y_control = mean(grandavg_control_12m.avg([1:12 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_control_12m.var([1:12 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_control_12m.var([1:12 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot adult data
figure;
time_x = grandavg_baseline_ad.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_baseline_ad.avg([4 7:8 10:12 14:18 24:25 27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_baseline_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_baseline_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("Adult ERP (Music vs Shuffled)");

% Average over all PMs for control condition
y_control = mean(grandavg_control_ad.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_control_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_control_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

%% pitch condititions
% Plot 3m data
figure;
time_x = grandavg_highvoice_3m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_highvoice_3m.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_highvoice_3m.var([1:13 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_highvoice_3m.var([1:13 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("3m ERP (High vs Low Pitch)");

% Average over all PMs for control condition
y_control = mean(grandavg_lowbass_3m.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_lowbass_3m.var([1:13 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_lowbass_3m.var([1:13 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot 6m data
figure;
time_x = grandavg_highvoice_6m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_highvoice_6m.avg([1:18 24:27 ],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_highvoice_6m.var([1:18 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_highvoice_6m.var([1:18 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("6m ERP (High vs Low Pitch)");

% Average over all PMs for control condition
y_control = mean(grandavg_lowbass_6m.avg([1:18 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_lowbass_6m.var([1:18 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_lowbass_6m.var([1:18 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot 12m data
figure;
time_x = grandavg_highvoice_12m.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_highvoice_12m.avg([1:12 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_highvoice_12m.var([1:12 24:27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_highvoice_12m.var([1:12 24:27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("12m ERP (High vs Low Pitch)");

% Average over all PMs for control condition
y_control = mean(grandavg_lowbass_12m.avg([1:12 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_lowbass_12m.var([1:12 24:27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_lowbass_12m.var([1:12 24:27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

% Plot adult data
figure;
time_x = grandavg_highvoice_ad.time(1,:);

% Average over all PMs for baseline condition
y_baseline = mean(grandavg_highvoice_ad.avg([4 7:8 10:12 14:18 24:25 27],:), 1); % Averaging over the 1st dimension (PMs)
y_baseline_upper = y_baseline + (mean(sqrt(grandavg_highvoice_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24)); 
y_baseline_lower = y_baseline - (mean(sqrt(grandavg_highvoice_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));

% Plot baseline
fill([time_x, fliplr(time_x)], [y_baseline_upper, fliplr(y_baseline_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(time_x, y_baseline, 'b');
title("Adult ERP (High vs Low Pitch)");

% Average over all PMs for control condition
y_control = mean(grandavg_lowbass_ad.avg([1:13 24:27],:), 1); % Averaging over the 1st dimension (PMs)
y_control_upper = y_control + (mean(sqrt(grandavg_lowbass_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));
y_control_lower = y_control - (mean(sqrt(grandavg_lowbass_ad.var([4 7:8 10:12 14:18 24:25 27],:)), 1) / sqrt(24));

% Plot control
fill([time_x, fliplr(time_x)], [y_control_upper, fliplr(y_control_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_x, y_control, 'r');
hold off;

%% define the parameters for the statistical comparison
cfg                  = [];
cfg.latency          = [-0.100 0.500];

cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                   % permutation distribution.
cfg.minnbchan        = 2;          % minimum number of neighborhood channels that is
                                   % required for a selected sample to be included
                                   % in the clustering algorithm (default=0).
cfg.neighbours       = 'layout_sing_neighb.mat'; % see below
cfg.tail             = 1;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 1;
cfg.alpha            = 0.05;      % alpha level of the permutation test
cfg.numrandomization = 500;        % number of draws from the permutation distribution

Nsub=26;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1; 
cfg.uvar             = 2; 

stat_ad_m = ft_timelockstatistics(cfg, avg_baseline_ad{:}, avg_control_ad{:}); 
stat_ad_f = ft_timelockstatistics(cfg, avg_lowbass_ad{:}, avg_highvoice_ad{:}); 

Nsub=26;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1; 
cfg.uvar             = 2; 

stat_6m_m = ft_timelockstatistics(cfg, avg_highvoice_6m{:}, avg_control_6m{:}); 
stat_6m_f = ft_timelockstatistics(cfg, avg_highvoice_6m{:}, avg_lowbass_6m{:}); 


Nsub=26;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1; 
cfg.uvar             = 2; 

stat_3m_m = ft_timelockstatistics(cfg, avg_highvoice_3m{:}, avg_control_3m{:}); 
stat_3m_f = ft_timelockstatistics(cfg, avg_highvoice_3m{:}, avg_lowbass_3m{:}); 


Nsub=27;
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1; 
cfg.uvar             = 2; 

stat_12m_m = ft_timelockstatistics(cfg, avg_highvoice_12m{:}, avg_control_12m{:}); 
stat_12m_f = ft_timelockstatistics(cfg, avg_highvoice_12m{:}, avg_lowbass_12m{:}); 



%% look at clusters

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_MvsC_3m    = ft_math(cfg, grandavg_baseline_3m,grandavg_control_3m);
GA_MvsC_6m    = ft_math(cfg, grandavg_baseline_6m,grandavg_control_6m);
GA_MvsC_12m    = ft_math(cfg, grandavg_baseline_12m,grandavg_control_12m);
GA_HVvsLB_3m    = ft_math(cfg, grandavg_highvoice_3m,grandavg_lowbass_3m);
GA_HVvsLB_6m    = ft_math(cfg, grandavg_highvoice_6m,grandavg_lowbass_6m);
GA_HVvsLB_12m    = ft_math(cfg, grandavg_highvoice_12m,grandavg_lowbass_12m);
GA_MvsC_ad    = ft_math(cfg, grandavg_baseline_voice_ad,grandavg_control_voice_ad);
GA_HVvsLB_ad    = ft_math(cfg, grandavg_highvoice_ad,grandavg_lowbass_ad);


% ft_singleplotER(cfg,GA_MvsC_3m,GA_MvsC_6m,GA_MvsC_12m)
% ft_singleplotER(cfg,GA_BvsLB_3m,GA_BvsHV_3m)
% ft_singleplotER(cfg,GA_BvsLB_6m,GA_BvsHV_6m)
% ft_singleplotER(cfg,GA_BvsLB_12m,GA_BvsHV_12m)

%%
stat=stat_12m;
GA_MvsC=GA_MvsC_12m;

figure;
% define parameters for plotting
timestep      = 0.02; %(in seconds)
sampling_rate = 1000;
sample_count  = length(stat.time);
j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

% get relevant values
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.05);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% get relevant values
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.05);
neg       = ismember(stat.negclusterslabelmat, neg_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(GA_MvsC.label, stat.label);

% plot
for k = 1:31;
   cfg.figure     = subplot(5,6,k);
   cfg.xlim       = [j(k) j(k+1)];
   cfg.zlim       = [-5e-14 5e-14];
   pos_int        = zeros(numel(GA_MvsC.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = 'layout_sing.mat';
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, GA_MvsC);
end

%% 

cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'layout_sing';33
% cfg.layout = 'CTF151_helmet.mat';
cfg.alpha     = 0.05;
% cfg.zlim      = [-4 4];
cfg.highlight = 'on';
cfg.highlightchannel = find(stat_6m.mask);
% cfg.comment   = 'no';
cfg.parameter = 'stat';
% ft_singleplotER(cfg, stat_3m);%, stat_6m, stat_12m)
ft_clusterplot(cfg, stat_6m);
title('Nonparametric: significant with cluster-based multiple comparison correction')

%% test condition differences



cfg = [];
cfg.channel = {'F3','F4','Fz','FC3','FC4','FCz','C3','C4','Cz'};
% cfg.latency = [0.000 0.444];%'all';
cfg.latency = [0.104 0.227];%positive cluster
cfg.latency = [0.307 0.325];%positive cluster 2
% cfg.latency = [-0.037 0.087];%negative cluster

cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_12m = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});
grandavg_control_12m = ft_timelockgrandaverage(cfg, avg_control_12m{:});

cfg = [];
cfg.channel = {'F3','F4','Fz','FC3','FC4','FCz','C3','C4','Cz', 'FP2','CP3','CP4','Pz','P4','P3'};
% cfg.latency = [0.000 0.444];%'all';
cfg.latency = [0.116 0.284];%positive cluster
% cfg.latency = [-0.077 0.123];%negative cluster
cfg.latency = [0.178 0.332];%positive cluster

cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_highvoice_6m = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_baseline_6m = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});
grandavg_control_6m = ft_timelockgrandaverage(cfg, avg_control_6m{:});

cfg = [];
cfg.channel = {'all'};
% cfg.latency = [0.000 0.444];%'all';
cfg.latency = [0.177 0.305];% positive cluster
% cfg.latency = [-0.037 0.106];%negative cluster

cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_highvoice_3m = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_baseline_3m = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});
grandavg_control_3m = ft_timelockgrandaverage(cfg, avg_control_3m{:});



cfg=[];
ft_singleplotER(cfg,grandavg_baseline_3m,grandavg_lowbass_3m,grandavg_highvoice_3m)
legend('Baseline 3m','Low Bass 3m','High Voice 3m')
ft_singleplotER(cfg,grandavg_highvoice_6m)
legend('Baseline 6m','Low Bass 6m','High Voice 6m')
ft_singleplotER(cfg,grandavg_baseline_12m,grandavg_lowbass_12m,grandavg_highvoice_12m)
legend('Baseline 12m','Low Bass 12m','High Voice 12m')

%% prep for stats
table_12m=[grandavg_baseline_12m.individual;grandavg_control_12m.individual;...
    grandavg_lowbass_12m.individual;grandavg_highvoice_12m.individual];
table_6m=[grandavg_baseline_6m.individual;grandavg_control_6m.individual;...
    grandavg_lowbass_6m.individual;grandavg_highvoice_6m.individual];
table_3m=[grandavg_baseline_3m.individual;grandavg_control_3m.individual;...
    grandavg_lowbass_3m.individual;grandavg_highvoice_3m.individual];

table_12m_peak = squeeze(max(table_12m,[], 3));
table_3m_peak = squeeze(max(table_3m,[], 3));
table_6m_peak = squeeze(max(table_6m,[], 3));

[table_12m_peak, idx] = max(table_12m, [], 3);
[row, ~] = ind2sub(size(table_12m, 1), idx);
table_12m_lat = grandavg_baseline_12m.time(row);

[table_6m_peak, idx] = max(table_6m, [], 3);
[row, ~] = ind2sub(size(table_6m, 1), idx);
table_6m_lat = grandavg_baseline_6m.time(row);

[table_3m_peak, idx] = max(table_3m, [], 3);
[row, ~] = ind2sub(size(table_3m, 1), idx);
table_3m_lat = grandavg_baseline_3m.time(row);

[table_6m_hv_peak, idx] = max(grandavg_highvoice_6m.individual, [], 3);
[row, ~] = ind2sub(size(grandavg_highvoice_6m, 1), idx);
table_6m_hv_lat = grandavg_highvoice_6m.time(row);


% table_3m_avg = nanmean(table_3m,3);
% table_6m_avg = nanmean(table_6m,3);
% table_12m_avg = nanmean(table_12m,3);

table_12m_auc = trapz(abs(table_12m),3);
table_6m_auc = trapz(abs(table_6m),3);
table_3m_auc = trapz(abs(table_3m),3);

table_peak = [table_12m_peak;table_3m_peak;table_6m_peak];
table_auc = [table_12m_auc;table_3m_auc;table_6m_auc];
table_lat = [table_12m_lat;table_3m_lat;table_6m_lat];

t = array2table(table_auc,'VariableNames',grandavg_baseline_3m.label);
writetable(t, "D:\MUSICOM_EEG\PROC\ERP_all_auc.xls");

tp = array2table(table_peak,'VariableNames',grandavg_baseline_3m.label);
writetable(tp, "D:\MUSICOM_EEG\PROC\ERP_all_peak.xls");

tl = array2table(table_lat,'VariableNames',grandavg_baseline_3m.label);
writetable(tl, "D:\MUSICOM_EEG\PROC\ERP_all_lat.xls");

%% run permutation analysis
design = '3x3';
levels = {{'1','2','3','1','2','3','1','2','3'}',{'1','1','1','2','2','2','3','3','3'}'};
channels = {'C3','Cz','C4','F3','Fz','F4','FC3','FCz','FC4' };
field_oi = 'individual';
data = {grandavg_baseline_3m,grandavg_lowbass_3m,grandavg_highvoice_3m,...
    grandavg_baseline_6m,grandavg_lowbass_6m,grandavg_highvoice_6m,...
    grandavg_baseline_12m,grandavg_lowbass_12m,grandavg_highvoice_12m};
[ data_out_F, data_out_P ] = Giac_rmANOVA( data, design, levels, channels, field_oi );
