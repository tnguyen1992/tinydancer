%% Step 4 Tiny Dancers: Group analyses for frequencies
clear all
clc

bigbin=5;
smallbin=3;
% load all frequency data
cd('D:\MUSICOM_EEG\PROC\Step3_frequencydata\3m');
FileNames_3m = {dir(['*_FREQ' '.mat']).name} ;
for k=1:1:length(FileNames_3m)
    load(FileNames_3m{k})
    avg_baseline_ssep_3m{k}     = Trinh_SSEP_removebacknoise(freq_baseline, bigbin,smallbin);
    avg_control_ssep_3m{k}      = Trinh_SSEP_removebacknoise(freq_control, bigbin,smallbin);
    avg_highvoice_ssep_3m{k}    = Trinh_SSEP_removebacknoise(freq_highvoice, bigbin,smallbin);
    avg_lowbass_ssep_3m{k}      = Trinh_SSEP_removebacknoise(freq_lowbass, bigbin,smallbin);

    freq_mean = nanmean(nanmean([freq_baseline.powspctrm,freq_control.powspctrm, ...
        freq_highvoice.powspctrm,freq_lowbass.powspctrm]));
    
    avg_baseline_3m{k}     = Trinh_freq_normalize(freq_baseline,freq_mean);
    avg_control_3m{k}      = Trinh_freq_normalize(freq_control,freq_mean);
    avg_highvoice_3m{k}    = Trinh_freq_normalize(freq_highvoice,freq_mean);
    avg_lowbass_3m{k}      = Trinh_freq_normalize(freq_lowbass,freq_mean);
end

cd('D:\MUSICOM_EEG\PROC\Step3_frequencydata\6m');
FileNames_6m = {dir(['*_FREQ' '.mat']).name} ;
for k=1:1:length(FileNames_6m)
    load(FileNames_6m{k})
    avg_baseline_ssep_6m{k}     = Trinh_SSEP_removebacknoise(freq_baseline, bigbin,smallbin);
    avg_control_ssep_6m{k}      = Trinh_SSEP_removebacknoise(freq_control, bigbin,smallbin);
    avg_highvoice_ssep_6m{k}    = Trinh_SSEP_removebacknoise(freq_highvoice, bigbin,smallbin);
    avg_lowbass_ssep_6m{k}      = Trinh_SSEP_removebacknoise(freq_lowbass, bigbin,smallbin);

    freq_mean = nanmean(nanmean([freq_baseline.powspctrm,freq_control.powspctrm, ...
        freq_highvoice.powspctrm,freq_lowbass.powspctrm]));

    avg_baseline_6m{k}     = Trinh_freq_normalize(freq_baseline,freq_mean);
    avg_control_6m{k}      = Trinh_freq_normalize(freq_control,freq_mean);
    avg_highvoice_6m{k}    = Trinh_freq_normalize(freq_highvoice,freq_mean);
    avg_lowbass_6m{k}      = Trinh_freq_normalize(freq_lowbass,freq_mean);
end

cd('D:\MUSICOM_EEG\PROC\Step3_frequencydata\12m');
FileNames_12m = {dir(['*_FREQ' '.mat']).name} ;
for k=1:1:length(FileNames_12m)
    load(FileNames_12m{k})
    avg_baseline_ssep_12m{k}     = Trinh_SSEP_removebacknoise(freq_baseline, bigbin,smallbin);
    avg_control_ssep_12m{k}      = Trinh_SSEP_removebacknoise(freq_control, bigbin,smallbin);
    avg_highvoice_ssep_12m{k}    = Trinh_SSEP_removebacknoise(freq_highvoice, bigbin,smallbin);
    avg_lowbass_ssep_12m{k}      = Trinh_SSEP_removebacknoise(freq_lowbass, bigbin,smallbin);

    freq_mean = nanmean(nanmean([freq_baseline.powspctrm,freq_control.powspctrm, ...
        freq_highvoice.powspctrm,freq_lowbass.powspctrm]));
    avg_baseline_12m{k}     = Trinh_freq_normalize(freq_baseline,freq_mean);
    avg_control_12m{k}      = Trinh_freq_normalize(freq_control,freq_mean);
    avg_highvoice_12m{k}    = Trinh_freq_normalize(freq_highvoice,freq_mean);
    avg_lowbass_12m{k}      = Trinh_freq_normalize(freq_lowbass,freq_mean);
end

cd('D:\MUSICOM_EEG\PROC\Step3_frequencydata\adults');
Filenames = {dir(['*_FREQ' '.mat']).name} ;
for k=1:1:length(Filenames)
    load(Filenames{k})
    avg_baseline_ssep{k}     = Trinh_SSEP_removebacknoise(freq_baseline, bigbin,smallbin);
    avg_control_ssep{k}      = Trinh_SSEP_removebacknoise(freq_control, bigbin,smallbin);
    avg_highvoice_ssep{k}    = Trinh_SSEP_removebacknoise(freq_highvoice, bigbin,smallbin);
    avg_lowbass_ssep{k}      = Trinh_SSEP_removebacknoise(freq_lowbass, bigbin,smallbin);

    freq_mean = nanmean(nanmean([freq_baseline.powspctrm,freq_control.powspctrm, ...
        freq_highvoice.powspctrm,freq_lowbass.powspctrm]));
    avg_baseline{k}     = Trinh_freq_normalize(freq_baseline,freq_mean);
    avg_control{k}      = Trinh_freq_normalize(freq_control,freq_mean);
    avg_highvoice{k}    = Trinh_freq_normalize(freq_highvoice,freq_mean);
    avg_lowbass{k}      = Trinh_freq_normalize(freq_lowbass,freq_mean);
end

% calculate grand average for each condition
cfg = [];
cfg.channel = {'Fz','F3','F4','FCz','FC3','FC4','C3','C4','Cz'};
% cfg.channel='all';
cfg.keepindividual = 'yes';
cfg.foilim         = [1.5 3];

grandavg_baseline_3m = ft_freqgrandaverage(cfg, avg_baseline_3m{:});
grandavg_control_3m = ft_freqgrandaverage(cfg, avg_control_3m{:});
grandavg_highvoice_3m = ft_freqgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_freqgrandaverage(cfg, avg_lowbass_3m{:});

grandavg_baseline_ssep_3m = ft_freqgrandaverage(cfg, avg_baseline_ssep_3m{:});
grandavg_control_ssep_3m = ft_freqgrandaverage(cfg, avg_control_ssep_3m{:});
grandavg_highvoice_ssep_3m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_3m{:});
grandavg_lowbass_ssep_3m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_3m{:});

grandavg_baseline_6m = ft_freqgrandaverage(cfg, avg_baseline_6m{:});
grandavg_control_6m = ft_freqgrandaverage(cfg, avg_control_6m{:});
grandavg_highvoice_6m = ft_freqgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_freqgrandaverage(cfg, avg_lowbass_6m{:});

grandavg_baseline_ssep_6m = ft_freqgrandaverage(cfg, avg_baseline_ssep_6m{:});
grandavg_control_ssep_6m = ft_freqgrandaverage(cfg, avg_control_ssep_6m{:});
grandavg_highvoice_ssep_6m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_6m{:});
grandavg_lowbass_ssep_6m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_6m{:});

grandavg_baseline_12m = ft_freqgrandaverage(cfg, avg_baseline_12m{:});
grandavg_control_12m = ft_freqgrandaverage(cfg, avg_control_12m{:});
grandavg_highvoice_12m = ft_freqgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_freqgrandaverage(cfg, avg_lowbass_12m{:});

grandavg_baseline_ssep_12m = ft_freqgrandaverage(cfg, avg_baseline_ssep_12m{:});
grandavg_control_ssep_12m = ft_freqgrandaverage(cfg, avg_control_ssep_12m{:});
grandavg_highvoice_ssep_12m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_12m{:});
grandavg_lowbass_ssep_12m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_12m{:});

grandavg_baseline = ft_freqgrandaverage(cfg, avg_baseline{:});
grandavg_control = ft_freqgrandaverage(cfg, avg_control{:});
grandavg_highvoice = ft_freqgrandaverage(cfg, avg_highvoice{:});
grandavg_lowbass = ft_freqgrandaverage(cfg, avg_lowbass{:});

grandavg_baseline_ssep = ft_freqgrandaverage(cfg, avg_baseline_ssep{:});
grandavg_control_ssep = ft_freqgrandaverage(cfg, avg_control_ssep{:});
grandavg_highvoice_ssep = ft_freqgrandaverage(cfg, avg_highvoice_ssep{:});
grandavg_lowbass_ssep = ft_freqgrandaverage(cfg, avg_lowbass_ssep{:});

%%
% for k=1:1:length(FileNames)
cfg = [];
cfg.channel      = {'Fz','F3','F4','FCz','FC3','FC4','C3','C4','Cz'}; % Fz','F3','F4','FCz','FC3','FC4',
cfg.channel      = {'all'};

cfg.layout       = 'layout_sing.mat';
cfg.xlim         = [1.5 3];
% cfg.ylim         = [-20000 40000];
% cfg.parameter    = 'powspctrm';
% ft_singleplotER(cfg, grandavg_baseline_ssep_3m,grandavg_control_ssep_3m);
% ft_singleplotER(cfg, grandavg_baseline_ssep_6m,grandavg_control_ssep_6m);
% ft_singleplotER(cfg, grandavg_baseline_ssep_12m,grandavg_control_ssep_12m);
% cfg.layout       = 'lay_biosemi64.mat';
% ft_singleplotER(cfg, grandavg_baseline_ssep,grandavg_control_ssep);
    
cfg.layout       = 'layout_sing.mat';
ft_singleplotER(cfg, grandavg_baseline_3m,grandavg_control_3m);
ft_singleplotER(cfg, grandavg_baseline_6m,grandavg_control_6m);
ft_singleplotER(cfg, grandavg_baseline_12m,grandavg_control_12m);
cfg.layout       = 'lay_biosemi64.mat';
ft_singleplotER(cfg, grandavg_baseline,grandavg_control);


%%
cfg.layout       = 'layout_sing.mat';
ft_singleplotER(cfg, grandavg_lowbass_ssep_3m,grandavg_highvoice_ssep_3m);
ft_singleplotER(cfg, grandavg_lowbass_ssep_6m,grandavg_highvoice_ssep_6m);
ft_singleplotER(cfg, grandavg_lowbass_ssep_12m,grandavg_highvoice_ssep_12m);
cfg.layout       = 'lay_biosemi64.mat';
ft_singleplotER(cfg, grandavg_lowbass_ssep,grandavg_highvoice_ssep);
%     savefig([num2str(k) '_12m_td_lola_plot.fig']);
%     close
% end

%%
% cfg.layout       = 'layout_sing.mat';
ft_singleplotER(cfg, grandavg_lowbass_3m,grandavg_highvoice_3m);
ft_singleplotER(cfg, grandavg_lowbass_6m,grandavg_highvoice_6m);
ft_singleplotER(cfg, grandavg_lowbass_12m,grandavg_highvoice_12m);
cfg.layout       = 'lay_biosemi64.mat';
ft_singleplotER(cfg, grandavg_lowbass,grandavg_highvoice);
%% music vs shuffle
freq=grandavg_baseline_ssep_3m.freq;
m3_ssep=[squeeze(mean(grandavg_baseline_ssep_3m.powspctrm,1:2)),squeeze(mean(grandavg_control_ssep_3m.powspctrm,1:2))];
n_baseline = size(grandavg_baseline_ssep_3m.powspctrm, 1); % Number of subjects
n_control = size(grandavg_control_ssep_3m.powspctrm, 1);
m3_ssep_se = [squeeze(std(mean(grandavg_baseline_ssep_3m.powspctrm,2),1)) / sqrt(n_baseline), ...
              squeeze(std(mean(grandavg_control_ssep_3m.powspctrm,2),1)) / sqrt(n_control)];
m6_ssep=[squeeze(mean(grandavg_baseline_ssep_6m.powspctrm,1:2)),squeeze(mean(grandavg_control_ssep_6m.powspctrm,1:2))];
n_baseline = size(grandavg_baseline_ssep_6m.powspctrm, 1); % Number of subjects
n_control = size(grandavg_control_ssep_6m.powspctrm, 1);
m6_ssep_se = [squeeze(std(mean(grandavg_baseline_ssep_6m.powspctrm,2),1)) / sqrt(n_baseline), ...
              squeeze(std(mean(grandavg_control_ssep_6m.powspctrm,2),1)) / sqrt(n_control)];
m12_ssep=[squeeze(mean(grandavg_baseline_ssep_12m.powspctrm,1:2)),...
    squeeze(mean(grandavg_control_ssep_12m.powspctrm,1:2))];
n_baseline = size(grandavg_baseline_ssep_12m.powspctrm, 1); % Number of subjects
n_control = size(grandavg_control_ssep_12m.powspctrm, 1);
m12_ssep_se = [squeeze(std(mean(grandavg_baseline_ssep_12m.powspctrm,2),1)) / sqrt(n_baseline), ...
              squeeze(std(mean(grandavg_control_ssep_12m.powspctrm,2),1)) / sqrt(n_control)];
mad_ssep=[squeeze(mean(grandavg_baseline_ssep.powspctrm,1:2)),...
    squeeze(mean(grandavg_control_ssep.powspctrm,1:2))];
n_baseline = size(grandavg_baseline_ssep.powspctrm, 1); % Number of subjects
n_control = size(grandavg_control_ssep.powspctrm, 1);
mad_ssep_se = [squeeze(std(mean(grandavg_baseline_ssep.powspctrm,2),1)) / sqrt(n_baseline), ...
              squeeze(std(mean(grandavg_control_ssep.powspctrm,2),1)) / sqrt(n_control)];



subplot(4,1,1) 
bar_handle = plot(freq, m3_ssep); 
hold on;

% Get bar positions for error bars
ngroups = size(m3_ssep, 2);
nbars = size(m3_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m3_ssep(:, i), m3_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-20000,60000]);

subplot(4,1,2) 
bar_handle = bar(freq, m6_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(m6_ssep, 2);
nbars = size(m6_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m6_ssep(:, i), m6_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-20000,80000]);

subplot(4,1,3) 
bar_handle = bar(freq, m12_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(m12_ssep, 2);
nbars = size(m12_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m12_ssep(:, i), m12_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-15000,50000]);

subplot(4,1,4) 
bar_handle = bar(freq, mad_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(mad_ssep, 2);
nbars = size(mad_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, mad_ssep(:, i), mad_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-4000,8000]);
%% pitch
figure;
freq=grandavg_highvoice_ssep_3m.freq;
m3_ssep=[squeeze(mean(grandavg_highvoice_ssep_3m.powspctrm,1:2)),squeeze(mean(grandavg_lowbass_ssep_3m.powspctrm,1:2))];
n_highvoice = size(grandavg_highvoice_ssep_3m.powspctrm, 1); % Number of subjects
n_lowbass = size(grandavg_lowbass_ssep_3m.powspctrm, 1);
m3_ssep_se = [squeeze(std(mean(grandavg_highvoice_ssep_3m.powspctrm,2),1)) / sqrt(n_highvoice), ...
              squeeze(std(mean(grandavg_lowbass_ssep_3m.powspctrm,2),1)) / sqrt(n_lowbass)];
m6_ssep=[squeeze(mean(grandavg_highvoice_ssep_6m.powspctrm,1:2)),squeeze(mean(grandavg_lowbass_ssep_6m.powspctrm,1:2))];
n_highvoice = size(grandavg_highvoice_ssep_6m.powspctrm, 1); % Number of subjects
n_lowbass = size(grandavg_lowbass_ssep_6m.powspctrm, 1);
m6_ssep_se = [squeeze(std(mean(grandavg_highvoice_ssep_6m.powspctrm,2),1)) / sqrt(n_highvoice), ...
              squeeze(std(mean(grandavg_lowbass_ssep_6m.powspctrm,2),1)) / sqrt(n_lowbass)];
m12_ssep=[squeeze(mean(grandavg_highvoice_ssep_12m.powspctrm,1:2)),...
    squeeze(mean(grandavg_lowbass_ssep_12m.powspctrm,1:2))];
n_highvoice = size(grandavg_highvoice_ssep_12m.powspctrm, 1); % Number of subjects
n_lowbass = size(grandavg_lowbass_ssep_12m.powspctrm, 1);
m12_ssep_se = [squeeze(std(mean(grandavg_highvoice_ssep_12m.powspctrm,2),1)) / sqrt(n_highvoice), ...
              squeeze(std(mean(grandavg_lowbass_ssep_12m.powspctrm,2),1)) / sqrt(n_lowbass)];
mad_ssep=[squeeze(mean(grandavg_highvoice_ssep.powspctrm,1:2)),...
    squeeze(mean(grandavg_lowbass_ssep.powspctrm,1:2))];
n_highvoice = size(grandavg_highvoice_ssep.powspctrm, 1); % Number of subjects
n_lowbass = size(grandavg_lowbass_ssep.powspctrm, 1);
mad_ssep_se = [squeeze(std(mean(grandavg_highvoice_ssep.powspctrm,2),1)) / sqrt(n_highvoice), ...
              squeeze(std(mean(grandavg_lowbass_ssep.powspctrm,2),1)) / sqrt(n_lowbass)];



subplot(4,1,1) 
bar_handle = bar(freq, m3_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(m3_ssep, 2);
nbars = size(m3_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m3_ssep(:, i), m3_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-20000,60000]);

subplot(4,1,2) 
bar_handle = bar(freq, m6_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(m6_ssep, 2);
nbars = size(m6_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m6_ssep(:, i), m6_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-20000,90000]);

subplot(4,1,3) 
bar_handle = bar(freq, m12_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(m12_ssep, 2);
nbars = size(m12_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, m12_ssep(:, i), m12_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-15000,60000]);

subplot(4,1,4) 
bar_handle = bar(freq, mad_ssep, 'grouped'); 
hold on;

% Get bar positions for error bars
ngroups = size(mad_ssep, 2);
nbars = size(mad_ssep, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5)); % Control group width

for i = 1:ngroups
    x = bar_handle(i).XData + (bar_handle(i).XOffset); % Bar x positions
    errorbar(x, mad_ssep(:, i), mad_ssep_se(:, i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

hold off;
ylim([-4000,10000]);

%% cluster based
cfg         = [];
cfg.frequency = [2.25];
cfg.channel          = {'Fz','F3','F4','FCz','FC3','FC4','C3','C4','Cz'};
cfg.avgoverchan      = 'yes';
cfg.method           = 'stats'; % montecarlo
cfg.statistic        = 'paired-ttest';
cfg.correctm         = 'none';
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.alpha            = 0.05;
% cfg.numrandomization = 1000;

Nsubj  = 26;
design = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
[stat_3m] = ft_freqstatistics(cfg, avg_baseline_ssep_3m{:},avg_control_ssep_3m{:});
[stat_6m] = ft_freqstatistics(cfg, avg_baseline_ssep_6m{:},avg_control_ssep_6m{:});
[stat_ad] = ft_freqstatistics(cfg, avg_baseline_ssep{:},avg_control_ssep{:});

[stat_6m_pitch] = ft_freqstatistics(cfg, avg_highvoice_ssep_6m{:},avg_lowbass_ssep_6m{:});

Nsubj  = 28;
design = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
[stat_12m] = ft_freqstatistics(cfg, avg_baseline_ssep_12m{:},avg_control_ssep_12m{:});
[stat_12m_pitch] = ft_freqstatistics(cfg, avg_lowbass_ssep_12m{:},avg_highvoice_ssep_12m{:});
%%
Nsubj  = 26;
design = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
[stat_3m] = ft_freqstatistics(cfg, avg_baseline_3m{:},avg_control_3m{:});
[stat_6m] = ft_freqstatistics(cfg, avg_baseline_6m{:},avg_control_6m{:});
[stat_ad] = ft_freqstatistics(cfg, avg_baseline{:},avg_control{:});

[stat_6m_pitch] = ft_freqstatistics(cfg, avg_highvoice_6m{:},avg_lowbass_6m{:});

Nsubj  = 28;
design = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
[stat_12m] = ft_freqstatistics(cfg, avg_baseline_12m{:},avg_control_12m{:});
[stat_12m_pitch] = ft_freqstatistics(cfg, avg_lowbass_12m{:},avg_highvoice_12m{:});


% bar(stat.freq,stat.prob)

% cfg = [];
% cfg.alpha  = 0.05;
% % cfg.parameter = 'freq';
% cfg.zlim   = [-1e-27 1e-27];
% cfg.layout = 'layout_sing.mat';
% ft_clusterplot(cfg, stat);
% % %% compare frequency between conditions
% % cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GA_BLvsCO    = ft_math(cfg, grandavg_baseline, grandavg_control);
%
% figure;
% % define parameters for plotting
% timestep      = 0.05; %(in seconds)
% sampling_rate = 1000;
% sample_count  = length(stat.time);
% j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
%
% % get relevant values
% pos_cluster_pvals = [stat.posclusters(:).prob];
% pos_clust = find(pos_cluster_pvals < 0.025);
% pos       = ismember(stat.posclusterslabelmat, pos_clust);
%
% % First ensure the channels to have the same order in the average and in the statistical output.
% % This might not be the case, because ft_math might shuffle the order
% [i1,i2] = match_str(GA_FICvsFC.label, stat.label);
%
% % plot
% for k = 1:20;
%    cfg.figure     = subplot(4,5,k);
%    cfg.xlim       = [j(k) j(k+1)];
%    cfg.zlim       = [-5e-14 5e-14];
%    pos_int        = zeros(numel(GA_FICvsFC.label),1);
%    pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
%    cfg.highlight  = 'on';
%    cfg.highlightchannel = find(pos_int);
%    cfg.comment    = 'xlim';
%    cfg.commentpos = 'title';
%    cfg.layout     = 'CTF151_helmet.mat';
%    cfg.figure     = 'gca';
%    ft_topoplotER(cfg, GA_FICvsFC);
% end


%% extract power information for stats
% calculate grand average for each condition
cfg = [];
cfg.channel = {'Fz','F3','F4','FCz','FC3','FC4','C3','C4','Cz'};
cfg.keepindividual = 'yes';
cfg.foilim         = [2.2 2.3];

grandavg_baseline_3m = ft_freqgrandaverage(cfg, avg_baseline_3m{:});
grandavg_control_3m = ft_freqgrandaverage(cfg, avg_control_3m{:});
grandavg_highvoice_3m = ft_freqgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_freqgrandaverage(cfg, avg_lowbass_3m{:});

grandavg_baseline_ssep_3m = ft_freqgrandaverage(cfg, avg_baseline_ssep_3m{:});
grandavg_control_ssep_3m = ft_freqgrandaverage(cfg, avg_control_ssep_3m{:});
grandavg_highvoice_ssep_3m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_3m{:});
grandavg_lowbass_ssep_3m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_3m{:});

grandavg_baseline_6m = ft_freqgrandaverage(cfg, avg_baseline_6m{:});
grandavg_control_6m = ft_freqgrandaverage(cfg, avg_control_6m{:});
grandavg_highvoice_6m = ft_freqgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_freqgrandaverage(cfg, avg_lowbass_6m{:});

grandavg_baseline_ssep_6m = ft_freqgrandaverage(cfg, avg_baseline_ssep_6m{:});
grandavg_control_ssep_6m = ft_freqgrandaverage(cfg, avg_control_ssep_6m{:});
grandavg_highvoice_ssep_6m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_6m{:});
grandavg_lowbass_ssep_6m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_6m{:});

grandavg_baseline_12m = ft_freqgrandaverage(cfg, avg_baseline_12m{:});
grandavg_control_12m = ft_freqgrandaverage(cfg, avg_control_12m{:});
grandavg_highvoice_12m = ft_freqgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_freqgrandaverage(cfg, avg_lowbass_12m{:});

grandavg_baseline_ssep_12m = ft_freqgrandaverage(cfg, avg_baseline_ssep_12m{:});
grandavg_control_ssep_12m = ft_freqgrandaverage(cfg, avg_control_ssep_12m{:});
grandavg_highvoice_ssep_12m = ft_freqgrandaverage(cfg, avg_highvoice_ssep_12m{:});
grandavg_lowbass_ssep_12m = ft_freqgrandaverage(cfg, avg_lowbass_ssep_12m{:});

grandavg_baseline = ft_freqgrandaverage(cfg, avg_baseline{:});
grandavg_control = ft_freqgrandaverage(cfg, avg_control{:});
grandavg_highvoice = ft_freqgrandaverage(cfg, avg_highvoice{:});
grandavg_lowbass = ft_freqgrandaverage(cfg, avg_lowbass{:});

grandavg_baseline_ssep = ft_freqgrandaverage(cfg, avg_baseline_ssep{:});
grandavg_control_ssep = ft_freqgrandaverage(cfg, avg_control_ssep{:});
grandavg_highvoice_ssep = ft_freqgrandaverage(cfg, avg_highvoice_ssep{:});
grandavg_lowbass_ssep = ft_freqgrandaverage(cfg, avg_lowbass_ssep{:});

%%
powerbaselinedf=nanmean([grandavg_baseline_ssep_3m.powspctrm; ...
    grandavg_baseline_ssep_6m.powspctrm; ...
    grandavg_baseline_ssep_12m.powspctrm; ...
    grandavg_baseline_ssep.powspctrm],2);
powercontroldf=nanmean([grandavg_control_ssep_3m.powspctrm; ...
    grandavg_control_ssep_6m.powspctrm; ...
    grandavg_control_ssep_12m.powspctrm; ...
    grandavg_control_ssep.powspctrm],2);
powerlowbassdf=nanmean([grandavg_lowbass_ssep_3m.powspctrm; ...
    grandavg_lowbass_ssep_6m.powspctrm; ...
    grandavg_lowbass_ssep_12m.powspctrm; ...
    grandavg_lowbass_ssep.powspctrm],2);
powerhighvoicedf=nanmean([grandavg_highvoice_ssep_3m.powspctrm; ...
    grandavg_highvoice_ssep_6m.powspctrm; ...
    grandavg_highvoice_ssep_12m.powspctrm; ...
    grandavg_highvoice_ssep.powspctrm],2);

b=regexp([FileNames_3m,FileNames_6m,FileNames_12m],'\d+(\.)?(\d+)?','match');
ID=str2double([b{:}])';

age3=repmat(3,26,1);
age6=repmat(6,26,1);
age12=repmat(12,28,1);

age=[age3;age6;age12];


condition=repmat('music',80,1);
power=powerbaselinedf;
T_music=table(ID,age, condition, power);
condition=repmat('control',80,1);
power=powercontroldf;
T_control=table(ID,age, condition, power);

condition=repmat('lowbass',80,1);
power=powerlowbassdf;
T_lowbass=table(ID,age, condition, power);
condition=repmat('highvoice',80,1);
power=powerhighvoicedf;
T_highvoice=table(ID,age, condition, power);


filename="power_stats_22535.xlsx";
sheetname = 'music';
writetable(T_music, filename, 'Sheet', sheetname);
sheetname = 'control';
writetable(T_control, filename, 'Sheet', sheetname);
sheetname = 'lowbass';
writetable(T_lowbass, filename, 'Sheet', sheetname);
sheetname = 'highvoice';
writetable(T_highvoice, filename, 'Sheet', sheetname);