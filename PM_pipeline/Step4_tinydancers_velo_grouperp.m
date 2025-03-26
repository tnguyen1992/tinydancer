%% Step 4 Tiny Dancers: Group analyses for ERP
% load all ERP data

%%
addpath 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\functions'
addpath 'D:\MUSICOM_EEG/onsets'
addpath('D:\MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
ft_defaults;
%% Load data
clear
clc
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc\3m\');
respdata_mus=[];
respdata_freq=[];

FileNames = {dir(['*' '.mat']).name} ;
for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_baseline_3m{k}     = ERP_baseline;
    avg_control_3m{k}      = ERP_control;
    respdata_mus = [respdata_mus; k, 3, numreps];
    clear ERP_baseline ERP_control
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc_freq\3m\');
FileNames = {dir(['*' '.mat']).name} ;
for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_3m{k}    = ERP_highvoice;
    avg_lowbass_3m{k}      = ERP_lowbass;
    respdata_freq = [respdata_freq; k, 3, numreps];

    clear ERP_highvoice ERP_lowbass 
end



% Initialize output array
names_3m = [];

% Loop over each element in the cell array
for i = 1:numel(FileNames)
    % Use regular expression to extract numerical substrings from string
    matches = regexp(FileNames{i}, '\d+(\.\d+)?', 'match');

    % Convert numerical substrings to double precision values
    values = str2double(matches);

    % Append values to output array
    names_3m = [names_3m, values];
end

cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc\6m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})

    avg_baseline_6m{k}     = ERP_baseline;
    avg_control_6m{k}      = ERP_control;
    respdata_mus = [respdata_mus; k+100, 6, numreps];
    
    clear ERP_baseline ERP_control
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc_freq\6m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_6m{k}    = ERP_highvoice;
    avg_lowbass_6m{k}      = ERP_lowbass;
    respdata_freq = [respdata_freq; k+100, 6, numreps];

    clear ERP_highvoice ERP_lowbass
end



% Initialize output array
names_6m = [];

% Loop over each element in the cell array
for i = 1:numel(FileNames)
    % Use regular expression to extract numerical substrings from string
    matches = regexp(FileNames{i}, '\d+(\.\d+)?', 'match');

    % Convert numerical substrings to double precision values
    values = str2double(matches);

    % Append values to output array
    names_6m = [names_6m, values];
end

cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc\12m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})

    avg_baseline_12m{k}     = ERP_baseline;
    avg_control_12m{k}      = ERP_control;
    respdata_mus = [respdata_mus; k+200, 12, numreps];

    clear ERP_baseline ERP_control 
end
cd('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc_freq\12m\');
FileNames = {dir(['*' '.mat']).name} ;

for k = 1:1:length(FileNames)
    load(FileNames{k})
    avg_highvoice_12m{k}    = ERP_highvoice;
    avg_lowbass_12m{k}      = ERP_lowbass;
    respdata_freq = [respdata_freq; k+200, 12, numreps];

    clear ERP_highvoice ERP_lowbass
end


% Initialize output array
names_12m = [];

% Loop over each element in the cell array
for i = 1:numel(FileNames)
    % Use regular expression to extract numerical substrings from string
    matches = regexp(FileNames{i}, '\d+(\.\d+)?', 'match');

    % Convert numerical substrings to double precision values
    values = str2double(matches);

    % Append values to output array
    names_12m = [names_12m, values];
end

% dlmwrite('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\IIT_Postdoc\WP4\MATLAB/respdata_freq.csv',respdata_freq)
% dlmwrite('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\IIT_Postdoc\WP4\MATLAB/respdata_mus.csv',respdata_mus)



%% plot ERPs
cfg = [];
cfg.layout = 'layout_sing.mat';
cfg.channel = {'all'};

cfg.xlim = [0 0.6];
% cfg.ylim = [35 75];
ft_singleplotER(cfg,avg_baseline_3m{:});
ft_singleplotER(cfg,avg_control_3m{:});%,{:})


ft_singleplotER(cfg,avg_baseline_6m{:});%,{:})
ft_singleplotER(cfg,avg_control_6m{:});%,{:})


ft_singleplotER(cfg,avg_baseline_12m{:});%,{:})
ft_singleplotER(cfg,avg_control_12m{:});%,{:})


%%
ft_singleplotER(cfg,avg_highvoice_3m{:});
ft_singleplotER(cfg,avg_highvoice_6m{:});
ft_singleplotER(cfg,avg_highvoice_12m{:});
ft_singleplotER(cfg,avg_lowbass_3m{:});
ft_singleplotER(cfg,avg_lowbass_6m{:});
ft_singleplotER(cfg,avg_lowbass_12m{:});

%% calculate grand average for each condition
cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.1 0.6];
cfg.parameter = 'avg';
cfg.nanmean        = 'no';
cfg.normalizevar   =  'N';
grandavg_control_3m = ft_timelockgrandaverage(cfg, avg_control_3m{:});
grandavg_control_6m = ft_timelockgrandaverage(cfg, avg_control_6m{:});
grandavg_control_12m = ft_timelockgrandaverage(cfg, avg_control_12m{:});

grandavg_baseline_12m = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});

grandavg_highvoice_6m = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_baseline_6m = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});


grandavg_highvoice_3m = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_baseline_3m = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});

%%
cfg = [];
cfg.layout = 'layout_sing.mat';
cfg.channel = {'all'};
cfg.xlim = [0 0.6];

ft_singleplotER(cfg,grandavg_baseline_3m,grandavg_control_3m)
ft_singleplotER(cfg,grandavg_baseline_6m,grandavg_control_6m)
ft_singleplotER(cfg,grandavg_baseline_12m,grandavg_control_12m)
%%

ft_singleplotER(cfg,grandavg_lowbass_3m, grandavg_highvoice_3m)
ft_singleplotER(cfg,grandavg_lowbass_6m, grandavg_highvoice_6m)
ft_singleplotER(cfg,grandavg_lowbass_12m, grandavg_highvoice_12m)
% ft_singleplotER(cfg,avg_highvoice_12m{:});%,avg_control_12m{:})

%% in individidual PMs in each age group
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


%% in individidual PMs in each age group pitch
figure;
for k = 1:10
    subplot(5,2,k)
    plot(grandavg_lowbass_3m.time(1,:),grandavg_lowbass_3m.avg(k,:))
    %     ylim([20 150]);
    title(sprintf("3m Movement to Tones in BP %d", k))
    hold on
    plot(grandavg_baseline_3m.time(1,:),grandavg_highvoice_3m.avg(k,:))
    %     ylim([20 150]);
    hold off
end

figure;
for k = 1:10
    subplot(5,2,k)
    plot(grandavg_lowbass_6m.time(1,:),grandavg_lowbass_6m.avg(k,:))
    %     ylim([20 150]);
    title(sprintf("6m Movement to Tones in BP %d", k))
    hold on
    plot(grandavg_baseline_6m.time(1,:),grandavg_highvoice_6m.avg(k,:))
    %     ylim([20 150]);
    hold off
end
figure;
for k = 1:10
    subplot(5,2,k)
    plot(grandavg_lowbass_12m.time(1,:),grandavg_lowbass_12m.avg(k,:))
    %     ylim([20 150]);
    title(sprintf("12m Movement to Tones in BP %d", k))
    hold on
    plot(grandavg_baseline_12m.time(1,:),grandavg_highvoice_12m.avg(k,:))
    %     ylim([20 150]);
    hold off
end

%% in individidual PMs in each age group
for k = 1:10
    subplot(5,2,k)
    plot(grandavg_baseline_3m.time(1,:),grandavg_baseline_3m.avg(k,:))
    title(sprintf("Movement to Music in PM %d", k))
    hold on
    plot(grandavg_baseline_6m.time(1,:),grandavg_baseline_6m.avg(k,:))
    plot(grandavg_baseline_12m.time(1,:),grandavg_baseline_12m.avg(k,:))
    hold off
end

figure;
for k = 1:10
    subplot(5,2,k)
    plot(grandavg_baseline_3m.time(1,:),grandavg_control_3m.avg(k,:))
    title(sprintf("Movement to Control in PM %d", k))
    hold on
    plot(grandavg_baseline_6m.time(1,:),grandavg_control_6m.avg(k,:))
    plot(grandavg_baseline_12m.time(1,:),grandavg_control_12m.avg(k,:))
    hold off
end
%% PMs superimposed
subplot(3,1,1)
for k = 1:10
    plot(grandavg_baseline_3m.time(1,:),grandavg_baseline_3m.avg(k,:),'b')
    title(sprintf("3m Movement to Music in PMs"))
    hold on
    plot(grandavg_baseline_3m.time(1,:),grandavg_control_3m.avg(k,:),'r')

end
hold off
subplot(3,1,2)
for k = 1:10
    plot(grandavg_baseline_6m.time(1,:),grandavg_baseline_6m.avg(k,:),'b')
    title(sprintf("6m Movement to Music in PMs"))
    hold on
    plot(grandavg_baseline_6m.time(1,:),grandavg_control_6m.avg(k,:),'r')
end
hold off
subplot(3,1,3)
for k = 1:10
    plot(grandavg_baseline_12m.time(1,:),grandavg_baseline_12m.avg(k,:),'b')
    title(sprintf("12m Movement to Music in PMs"))
    hold on
    plot(grandavg_baseline_12m.time(1,:),grandavg_control_12m.avg(k,:),'r')
end
hold off


% figure
% subplot(3,1,1)
% for k = 1:10
%     plot(grandavg_baseline_3m.time(1,:),grandavg_baseline_3m_dm.avg(k,:),'b')
%     title(sprintf("3m Movement to Music in PMs"))
%     hold on
%     plot(grandavg_baseline_3m.time(1,:),grandavg_control_3m_dm.avg(k,:),'r')
% 
% end
% hold off
% subplot(3,1,2)
% for k = 1:10
%     plot(grandavg_baseline_6m.time(1,:),grandavg_baseline_6m_dm.avg(k,:),'b')
%     title(sprintf("6m Movement to Music in PMs"))
%     hold on
%     plot(grandavg_baseline_6m.time(1,:),grandavg_control_6m_dm.avg(k,:),'r')
% end
% hold off
% subplot(3,1,3)
% for k = 1:10
%     plot(grandavg_baseline_12m.time(1,:),grandavg_baseline_12m_dm.avg(k,:),'b')
%     title(sprintf("12m Movement to Music in PMs"))
%     hold on
%     plot(grandavg_baseline_12m.time(1,:),grandavg_control_12m_dm.avg(k,:),'r')
% end
% hold off

%% define the parameters for the statistical comparison
cfg                  = [];
cfg.latency          = [0.00 0.500];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha            = 0.05;      % alpha level of the permutation test
% cfg.channel      = {'PM1', 'PM6'};
cfg.avgoverchan      = 'yes';
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


Nsub=length(avg_baseline_3m);
cfg.design(1, 1:2*Nsub) = [ones(1,Nsub), ones(1,Nsub)*2]; % design matrix
cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.ivar             = 1;
cfg.uvar             = 2;

stat_3m_mus = ft_timelockstatistics(cfg, avg_baseline_3m{:}, avg_control_3m{:});
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




%%
plotPM=10;
plotavg=1;
% subplot(3,1,1)
for k = 1:plotPM
    indices_below_005 = stat_3m.prob(k,:) < 0.05;
    subplot(5,2,k)
    plot(grandavg_baseline_3m.time(1,14:29),grandavg_baseline_3m.avg(k,14:29))
    title(sprintf("3m Movement to Tones in PM %d", k))
    hold on
    plot(grandavg_baseline_3m.time(1,14:29), grandavg_baseline_3m.avg(k,14:29).*indices_below_005, 'r.', 'markersize', 10); % You can customize the marker style/color
    plot(grandavg_baseline_3m.time(1,14:29),grandavg_control_3m.avg(k,14:29))
    hold off
end

figure
% subplot(3,1,2)
for k = 1:plotPM
    indices_below_005 = stat_6m.prob(k,:) < 0.05;
    subplot(5,2,k)
    plot(grandavg_baseline_6m.time(1,14:29),grandavg_baseline_6m.avg(k,14:29))
    title(sprintf("6m Movement to Tones in PM %d", k))
    hold on
    plot(grandavg_baseline_6m.time(1,14:29), grandavg_baseline_6m.avg(k,14:29).*indices_below_005, 'r.', 'markersize', 10); % You can customize the marker style/color
    plot(grandavg_baseline_6m.time(1,14:29),grandavg_control_6m.avg(k,14:29))
    hold off
end

% subplot(3,1,3)
figure
for k = 1:plotPM
    indices_below_005 = stat_12m.prob(k,:) < 0.05;
    subplot(5,2,k)
    plot(grandavg_baseline_12m.time(1,14:29),grandavg_baseline_12m.avg(k,14:29))
    title(sprintf("12m Movement to Tones in PM %d", k))
    hold on
    plot(grandavg_baseline_12m.time(1,14:29), grandavg_baseline_12m.avg(k,14:29).*indices_below_005, 'r.', 'markersize', 10); % You can customize the marker style/color
    plot(grandavg_baseline_12m.time(1,14:29),grandavg_control_12m.avg(k,14:29))
    hold off
end

%% look at clusters

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_MvsC_3m    = ft_math(cfg, grandavg_music_3m,grandavg_control_3m);
GA_MvsC_6m    = ft_math(cfg, grandavg_music_6m,grandavg_control_6m);
GA_MvsC_12m    = ft_math(cfg, grandavg_music_12m,grandavg_control_12m);


stat=stat_3m;
GA_MvsC=GA_MvsC_3m;

figure;
% define parameters for plotting
timestep      = 0.04; %(in seconds)
sampling_rate = 25;
sample_count  = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
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
for k = 1:10;
    cfg.figure     = subplot(5,5,k);
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
cfg.layout    = 'layout_sing';
% cfg.layout = 'CTF151_helmet.mat';
cfg.alpha     = 0.05;
% cfg.zlim      = [-4 4];
cfg.highlight = 'on';
cfg.highlightchannel = find(stat_3m.mask);
% cfg.comment   = 'no';
cfg.parameter = 'stat';
ft_singleplotER(cfg, stat_3m);%, stat_6m, stat_12m)
ft_singleplotER(cfg, stat_6m);%, stat_6m, stat_12m)
ft_singleplotER(cfg, stat_12m);%, stat_6m, stat_12m)

title('Nonparametric: significant with cluster-based multiple comparison correction')

%% test condition differences
% on beat (40%)

cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.088 0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_12m_onb = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m_onb = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m_onb = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});
grandavg_control_12m_onb = ft_timelockgrandaverage(cfg, avg_control_12m{:});

cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.088 0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_highvoice_6m_onb = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m_onb = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_baseline_6m_onb = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});
grandavg_control_6m_onb = ft_timelockgrandaverage(cfg, avg_control_6m{:});

cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.088 0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_highvoice_3m_onb = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m_onb = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_baseline_3m_onb = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});
grandavg_control_3m_onb = ft_timelockgrandaverage(cfg, avg_control_3m{:});


% off beat 60% (88.8-222 ms)
% 12m
% get the first 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.222 -0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_12m_offb1 = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m_offb1 = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m_offb1 = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});
grandavg_control_12m_offb1 = ft_timelockgrandaverage(cfg, avg_control_12m{:});

% get the other 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [0.088 0.222];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_12m_offb = ft_timelockgrandaverage(cfg, avg_baseline_12m{:});
grandavg_highvoice_12m_offb = ft_timelockgrandaverage(cfg, avg_highvoice_12m{:});
grandavg_lowbass_12m_offb = ft_timelockgrandaverage(cfg, avg_lowbass_12m{:});
grandavg_control_12m_offb = ft_timelockgrandaverage(cfg, avg_control_12m{:});

% put both back together
grandavg_baseline_12m_offb.individual = squeeze(nanmean(cat(4,grandavg_baseline_12m_offb1.individual,grandavg_baseline_12m_offb.individual),4));
grandavg_highvoice_12m_offb.individual = squeeze(nanmean(cat(4,grandavg_highvoice_12m_offb1.individual,grandavg_highvoice_12m_offb.individual),4));
grandavg_lowbass_12m_offb.individual = squeeze(nanmean(cat(4,grandavg_lowbass_12m_offb1.individual,grandavg_lowbass_12m_offb.individual),4));
grandavg_control_12m_offb.individual = squeeze(nanmean(cat(4,grandavg_control_12m_offb1.individual,grandavg_control_12m_offb.individual),4));

% 6m
% get the first 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.222 -0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_6m_offb1 = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});
grandavg_highvoice_6m_offb1 = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m_offb1 = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_control_6m_offb1 = ft_timelockgrandaverage(cfg, avg_control_6m{:});

% get the other 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [0.088 0.222];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_6m_offb = ft_timelockgrandaverage(cfg, avg_baseline_6m{:});
grandavg_highvoice_6m_offb = ft_timelockgrandaverage(cfg, avg_highvoice_6m{:});
grandavg_lowbass_6m_offb = ft_timelockgrandaverage(cfg, avg_lowbass_6m{:});
grandavg_control_6m_offb = ft_timelockgrandaverage(cfg, avg_control_6m{:});

% put both back together
grandavg_baseline_6m_offb.individual = squeeze(nanmean(cat(4,grandavg_baseline_6m_offb1.individual,grandavg_baseline_6m_offb.individual),4));
grandavg_highvoice_6m_offb.individual = squeeze(nanmean(cat(4,grandavg_highvoice_6m_offb1.individual,grandavg_highvoice_6m_offb.individual),4));
grandavg_lowbass_6m_offb.individual = squeeze(nanmean(cat(4,grandavg_lowbass_6m_offb1.individual,grandavg_lowbass_6m_offb.individual),4));
grandavg_control_6m_offb.individual = squeeze(nanmean(cat(4,grandavg_control_6m_offb1.individual,grandavg_control_6m_offb.individual),4));

% 3m
% get the first 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [-0.222 -0.088];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_3m_offb1 = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});
grandavg_highvoice_3m_offb1 = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m_offb1 = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_control_3m_offb1 = ft_timelockgrandaverage(cfg, avg_control_3m{:});

% get the other 30 %
cfg = [];
cfg.channel = {'all'};
cfg.latency = [0.088 0.222];%'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
grandavg_baseline_3m_offb = ft_timelockgrandaverage(cfg, avg_baseline_3m{:});
grandavg_highvoice_3m_offb = ft_timelockgrandaverage(cfg, avg_highvoice_3m{:});
grandavg_lowbass_3m_offb = ft_timelockgrandaverage(cfg, avg_lowbass_3m{:});
grandavg_control_3m_offb = ft_timelockgrandaverage(cfg, avg_control_3m{:});

% put both back together
grandavg_baseline_3m_offb.individual = squeeze(nanmean(cat(4,grandavg_baseline_3m_offb1.individual,grandavg_baseline_3m_offb.individual),4));
grandavg_highvoice_3m_offb.individual = squeeze(nanmean(cat(4,grandavg_highvoice_3m_offb1.individual,grandavg_highvoice_3m_offb.individual),4));
grandavg_lowbass_3m_offb.individual = squeeze(nanmean(cat(4,grandavg_lowbass_3m_offb1.individual,grandavg_lowbass_3m_offb.individual),4));
grandavg_control_3m_offb.individual = squeeze(nanmean(cat(4,grandavg_control_3m_offb1.individual,grandavg_control_3m_offb.individual),4));

%% calculate the beat contrast

grandavg_baseline_12m_velocontrast = grandavg_baseline_12m_onb;
grandavg_baseline_12m_velocontrast.individual = (grandavg_baseline_12m_onb.individual-grandavg_baseline_12m_offb.individual)./(grandavg_baseline_12m_onb.individual+grandavg_baseline_12m_offb.individual);
grandavg_highvoice_12m_velocontrast = grandavg_highvoice_12m_onb;
grandavg_highvoice_12m_velocontrast.individual = (grandavg_highvoice_12m_onb.individual-grandavg_highvoice_12m_offb.individual)./(grandavg_highvoice_12m_onb.individual+grandavg_highvoice_12m_offb.individual);
grandavg_lowbass_12m_velocontrast = grandavg_lowbass_12m_onb;
grandavg_lowbass_12m_velocontrast.individual = (grandavg_lowbass_12m_onb.individual-grandavg_lowbass_12m_offb.individual)./(grandavg_lowbass_12m_onb.individual+grandavg_lowbass_12m_offb.individual);
grandavg_control_12m_velocontrast = grandavg_control_12m_onb;
grandavg_control_12m_velocontrast.individual = (grandavg_control_12m_onb.individual-grandavg_control_12m_offb.individual)./(grandavg_control_12m_onb.individual+grandavg_control_12m_offb.individual);

grandavg_baseline_6m_velocontrast = grandavg_baseline_6m_onb;
grandavg_baseline_6m_velocontrast.individual = (grandavg_baseline_6m_onb.individual-grandavg_baseline_6m_offb.individual)./(grandavg_baseline_6m_onb.individual+grandavg_baseline_6m_offb.individual);
grandavg_highvoice_6m_velocontrast = grandavg_highvoice_6m_onb;
grandavg_highvoice_6m_velocontrast.individual = (grandavg_highvoice_6m_onb.individual-grandavg_highvoice_6m_offb.individual)./(grandavg_highvoice_6m_onb.individual+grandavg_highvoice_6m_offb.individual);
grandavg_lowbass_6m_velocontrast = grandavg_lowbass_6m_onb;
grandavg_lowbass_6m_velocontrast.individual = (grandavg_lowbass_6m_onb.individual-grandavg_lowbass_6m_offb.individual)./(grandavg_lowbass_6m_onb.individual+grandavg_lowbass_6m_offb.individual);
grandavg_control_6m_velocontrast = grandavg_control_6m_onb;
grandavg_control_6m_velocontrast.individual = (grandavg_control_6m_onb.individual-grandavg_control_6m_offb.individual)./(grandavg_control_6m_onb.individual+grandavg_control_6m_offb.individual);

grandavg_baseline_3m_velocontrast = grandavg_baseline_3m_onb;
grandavg_baseline_3m_velocontrast.individual = (grandavg_baseline_3m_onb.individual-grandavg_baseline_3m_offb.individual)./(grandavg_baseline_3m_onb.individual+grandavg_baseline_3m_offb.individual);
grandavg_highvoice_3m_velocontrast = grandavg_highvoice_3m_onb;
grandavg_highvoice_3m_velocontrast.individual = (grandavg_highvoice_3m_onb.individual-grandavg_highvoice_3m_offb.individual)./(grandavg_highvoice_3m_onb.individual+grandavg_highvoice_3m_offb.individual);
grandavg_lowbass_3m_velocontrast = grandavg_lowbass_3m_onb;
grandavg_lowbass_3m_velocontrast.individual = (grandavg_lowbass_3m_onb.individual-grandavg_lowbass_3m_offb.individual)./(grandavg_lowbass_3m_onb.individual+grandavg_lowbass_3m_offb.individual);
grandavg_control_3m_velocontrast = grandavg_control_3m_onb;
grandavg_control_3m_velocontrast.individual = (grandavg_control_3m_onb.individual-grandavg_control_3m_offb.individual)./(grandavg_control_3m_onb.individual+grandavg_control_3m_offb.individual);





%% plotting
bodyparts={'head', 'torso', 'leftelbow','lefthand', 'rightelbow','righthand','leftknee', 'leftfoot', 'rightknee', 'rightfoot'};
tiledlayout(2,5);

for b=1:length(bodyparts)

    if strcmp(bodyparts{b},'head')
        columns_x=[1,3,5];
        columns_y=[2,4,6];
    elseif strcmp(bodyparts{b},'torso')
        columns_x=[7,9,11];
        columns_y=[8,10,12];
    elseif strcmp(bodyparts{b},'rightelbow')
        columns_x=[13];
        columns_y=[14];
    elseif strcmp(bodyparts{b},'righthand')
        columns_x=[15,17];
        columns_y=[16,18];
    elseif strcmp(bodyparts{b},'leftelbow')
        columns_x=[19];
        columns_y=[20];
    elseif strcmp(bodyparts{b},'lefthand')
        columns_x=[21,23];
        columns_y=[22,24];
    elseif strcmp(bodyparts{b},'rightknee')
        columns_x=[25];
        columns_y=[26];
    elseif strcmp(bodyparts{b},'rightfoot')
        columns_x=[27,29];
        columns_y=[28,30];
    elseif strcmp(bodyparts{b},'leftknee')
        columns_x=[31];
        columns_y=[32];
    elseif strcmp(bodyparts{b},'leftfoot')
        columns_x=[33,35];
        columns_y=[34,36];
    end

    % baseline
    % 12 m
    data=grandavg_baseline_12m;
    time = data.time;
    bodypart_12_music_x = squeeze(nanmean(data.avg(columns_x,:),1));
    bodypart_12_music_y = squeeze(nanmean(data.avg(columns_y,:),1));

    % 6 m
    data=grandavg_baseline_6m;
    time = data.time;
    bodypart_6_music_x = squeeze(nanmean(data.avg(columns_x,:),1));
    bodypart_6_music_y = squeeze(nanmean(data.avg(columns_y,:),1));

    % 3 m
    data=grandavg_baseline_3m;
    time = data.time;
    bodypart_3_music_x = squeeze(nanmean(data.avg(columns_x,:),1));
    bodypart_3_music_y = squeeze(nanmean(data.avg(columns_y,:),1));

    % put everything together

    nexttile;
    plot(time,bodypart_12_music_x,'b-')
    hold on
    plot(time,bodypart_12_music_y,'b--')
    plot(time,bodypart_6_music_x,'r-')
    plot(time,bodypart_6_music_y,'r--')
    plot(time,bodypart_3_music_x,'g-')
    plot(time,bodypart_3_music_y,'g--')

    xline(0, '-'); % Draw a red dashed line at x = 3
    xline(0.088,'--');
    xline(-0.088,'--');
    xline(-0.222,'--');
    xline(0.222,'--');
    xline(0.444, '-'); % Draw a red dashed line at x = 3
    xline(-0.444, '-'); % Draw a red dashed line at x = 3
    title(bodyparts{b})
end

%%
cfg=[];
cfg.channel={'all'};
% cfg.ylim=[0 0.2];
% ft_singleplotER(cfg,grandavg_baseline_3m_onb,grandavg_control_3m_onb,grandavg_lowbass_3m_onb,grandavg_highvoice_3m_onb)
% legend('Baseline 3m','Control 3m','Low Bass 3m','High Voice 3m')
% ft_singleplotER(cfg,grandavg_baseline_6m_onb,grandavg_control_6m_onb,grandavg_lowbass_6m_onb,grandavg_highvoice_6m_onb)
% legend('Baseline 6m','Control 6m','Low Bass 6m','High Voice 6m')
% ft_singleplotER(cfg,grandavg_baseline_12m_onb,grandavg_control_12m_onb,grandavg_lowbass_12m_onb,grandavg_highvoice_12m_onb)
% legend('Baseline 12m','Control 12m','Low Bass 12m','High Voice 12m')
%
% ft_singleplotER(cfg,grandavg_baseline_3m_offb,grandavg_control_3m_offb,grandavg_lowbass_3m_offb,grandavg_highvoice_3m_offb)
% legend('Baseline 3m','Control 3m','Low Bass 3m','High Voice 3m')
% ft_singleplotER(cfg,grandavg_baseline_6m_offb,grandavg_control_6m_offb,grandavg_lowbass_6m_offb,grandavg_highvoice_6m_offb)
% legend('Baseline 6m','Control 6m','Low Bass 6m','High Voice 6m')
% ft_singleplotER(cfg,grandavg_baseline_12m_offb,grandavg_control_12m_offb,grandavg_lowbass_12m_offb,grandavg_highvoice_12m_offb)
% legend('Baseline 12m','Control 12m','Low Bass 12m','High Voice 12m')

ft_singleplotER(cfg,grandavg_baseline_3m_velocontrast,grandavg_control_3m_velocontrast,grandavg_lowbass_3m_velocontrast,grandavg_highvoice_3m_velocontrast)
legend('Baseline 3m','Control 3m','Low Bass 3m','High Voice 3m')
ft_singleplotER(cfg,grandavg_baseline_6m_velocontrast,grandavg_control_6m_velocontrast,grandavg_lowbass_6m_velocontrast,grandavg_highvoice_6m_velocontrast)
legend('Baseline 6m','Control 6m','Low Bass 6m','High Voice 6m')
ft_singleplotER(cfg,grandavg_baseline_12m_velocontrast,grandavg_control_12m_velocontrast,grandavg_lowbass_12m_velocontrast,grandavg_highvoice_12m_velocontrast)
legend('Baseline 12m','Control 12m','Low Bass 12m','High Voice 12m')

%% prep for stats

table_12m=[grandavg_baseline_12m_velocontrast.individual;grandavg_lowbass_12m_velocontrast.individual;grandavg_highvoice_12m_velocontrast.individual;grandavg_control_12m_velocontrast.individual];
table_3m=[grandavg_baseline_6m_velocontrast.individual;grandavg_lowbass_6m_velocontrast.individual;grandavg_highvoice_6m_velocontrast.individual;grandavg_control_6m_velocontrast.individual];
table_6m=[grandavg_baseline_3m_velocontrast.individual;grandavg_lowbass_3m_velocontrast.individual;grandavg_highvoice_3m_velocontrast.individual;grandavg_control_3m_velocontrast.individual];

table_12m_avg = nanmean(table_12m,3);
table_3m_avg = nanmean(table_3m,3);
table_6m_avg = nanmean(table_6m,3);

% insert the names
names_12m=repmat(names_12m,1,4).';
names_6m=repmat(names_6m,1,4).';
names_3m=repmat(names_3m,1,4).';

names=[names_12m;names_3m;names_6m];


% insert the conditions
condition={'baseline','highvoice','lowbass','control'};
% Initialize a new cell array to hold the repeated elements
repeated_conditions_12m = cell(1, numel(condition) * 28);

% Loop over each element in the condition cell array
for i = 1:numel(condition)
    % Repeat the current element 28 times using repmat
    repeated_element = repmat(condition(i), 1, 28);

    % Insert the repeated element into the new cell array
    start_index = (i-1) * 28 + 1;
    end_index = i * 28;
    repeated_conditions_12m(start_index:end_index) = repeated_element;
end

% Initialize a new cell array to hold the repeated elements
repeated_conditions_3m = cell(1, numel(condition) * 28);

% Loop over each element in the condition cell array
for i = 1:numel(condition)
    % Repeat the current element 28 times using repmat
    repeated_element = repmat(condition(i), 1, 28);

    % Insert the repeated element into the new cell array
    start_index = (i-1) * 28 + 1;
    end_index = i * 28;
    repeated_conditions_3m(start_index:end_index) = repeated_element;
end


% Initialize a new cell array to hold the repeated elements
repeated_conditions_6m = cell(1, numel(condition) * 25);

% Loop over each element in the condition cell array
for i = 1:numel(condition)
    % Repeat the current element 25 times using repmat
    repeated_element = repmat(condition(i), 1, 25);

    % Insert the repeated element into the new cell array
    start_index = (i-1) * 25 + 1;
    end_index = i * 25;
    repeated_conditions_6m(start_index:end_index) = repeated_element;
end

repeated_conditions = [repeated_conditions_12m.';repeated_conditions_3m.';repeated_conditions_6m.'];
% Convert repeated_conditions to a categorical array
repeated_conditions = categorical(repeated_conditions);

table_all = [table_12m_avg;table_3m_avg;table_6m_avg];
t1 = table(names, repeated_conditions,'VariableNames',{'ID','condition'});
t2 = array2table(table_all,'VariableNames', ERP_music.label);
%     {'lefteye_x', 'lefteye_y', 'righteye_x', 'righteye_y', 'mouth_x', 'mouth_y', 'leftshoulder_x', 'leftshoulder_y', 'rightshoulder_x', 'rightshoulder_y', 'chest_x', 'chest_y', 'leftelbow_x', 'leftelbow_y', 'leftwrist_x', 'leftwrist_y', 'lefthand_x', 'lefthand_y', 'rightelbow_x', 'rightelbow_y', 'rightwrist_x', 'rightwrist_y', 'righthand_x', 'righthand_y', 'leftknee_x', 'leftknee_y', 'leftankle_x', 'leftankle_y', 'leftfoot_x', 'leftfoot_y', 'rightknee_x', 'rightknee_y', 'rightankle_x', 'rightankle_y', 'rightfoot_x', 'rightfoot_y'});

t=[t1 t2];

writetable(t, "D:\TinyDancers_Movement\fieldtrip_velo\PROC\movement_beatcontrast_output.xls");
%% run stat analysis
design = '3x3';
levels = {{'1','2','3','1','2','3','1','2','3'}',{'1','1','1','2','2','2','3','3','3'}'};
channels = {'lefthand','righthand','leftfoot','rightfoot','leftwrist','rightwrist','leftankle','rightankle' };
field_oi = 'individual';
data = {grandavg_baseline_3m,grandavg_lowbass_3m,grandavg_highvoice_3m,...
    grandavg_baseline_6m,grandavg_lowbass_6m,grandavg_highvoice_6m,...
    grandavg_baseline_12m,grandavg_lowbass_12m,grandavg_highvoice_12m};
[ data_out_F, data_out_P ] = Giac_rmANOVA( data, design, levels, channels, field_oi );
