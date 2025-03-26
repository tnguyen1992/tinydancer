%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear 
clc

addpath 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\functions'
addpath 'D:\MUSICOM_EEG/onsets'
addpath('D:\MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
ft_defaults;
%% load the amplitude envelopes
load('D:\MUSICOM_EEG\SideProject-idyom-trf\Stimuli\StimAcoustics_movement.mat')

hopp_music=env{1,1};
hopp_control=env{3,1};
hopp_lp=env{2,1};
hopp_hp=env{5,1};

lola_music=env{6,1};
lola_control=env{8,1};
lola_lp=env{7,1};
lola_hp=env{10,1};

hopp_music_time = 0:0.04:21.32;
[pks_hb,locs_hb] = findpeaks(hopp_music,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_hb);
% indices_above_median = pks_hb > median_value;
% locs_hb = locs_hb(indices_above_median);

% plot(hopp_music_time,hopp_music,hopp_music_time(locs_hb),pks_hb,'or')
[pks_lb,locs_lb] = findpeaks(lola_music,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% plot(hopp_music_time,lola_music,hopp_music_time(locs_lb),pks_lb,'or')
% median_value = median(pks_lb);
% indices_above_median = pks_lb > median_value;
% locs_lb = locs_lb(indices_above_median);


[pks_hc,locs_hc] = findpeaks(hopp_control,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_hc);
% indices_above_median = pks_hc > median_value;
% locs_hc = locs_hc(indices_above_median);

% plot(hopp_music_time,hopp_control,hopp_music_time(locs_hc),pks_hc,'or')
[pks_lc,locs_lc] = findpeaks(lola_control,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_lc);
% indices_above_median = pks_lc > median_value;
% locs_lc = locs_lc(indices_above_median);
% plot(hopp_music_time,lola_control,hopp_music_time(locs_lc),pks_lc,'or')

[pks_hh,locs_hh] = findpeaks(hopp_hp,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_hh);
% indices_above_median = pks_hh > median_value;
% locs_hh = locs_hh(indices_above_median);

[pks_lh,locs_lh] = findpeaks(lola_hp,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_lh);
% indices_above_median = pks_lh > median_value;
% locs_lh = locs_lh(indices_above_median);

[pks_hl,locs_hl] = findpeaks(hopp_lp,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_hl);
% indices_above_median = pks_hl > median_value;
% locs_hl = locs_hl(indices_above_median);

[pks_ll,locs_ll] = findpeaks(lola_lp,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% median_value = median(pks_ll);
% indices_above_median = pks_ll > median_value;
% locs_ll = locs_ll(indices_above_median);


%% access preprocessed files
group={'3m','6m','12m'};
for g=1:length(group)
    cd(['D:\TinyDancers_Movement\fieldtrip_velo\PM\' group{g}]);
    FileNames = {dir(['*' '.mat']).name} ;

    for k=1:1:length(FileNames)

        % load data
        load(FileNames{k});
        baby_name       = erase(FileNames{k},'.mat');

        cfg=[];        
        data_baseline = ft_preprocessing(cfg,data_baseline);
        data_control = ft_preprocessing(cfg,data_control);
        data_highvoice = ft_preprocessing(cfg,data_highvoice);
        data_lowbass = ft_preprocessing(cfg,data_lowbass);

      


        %% Epoching

        % Reference condition
        if length(unique(data_baseline.trialinfo))>1
            [data_hbaseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_hb,3,[0.48 3.48]);
            [data_lbaseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_lb,7,[0.48 3.48]);
            cfg                 = [];
            cfg.keepsampleinfo  ='no';
            data_baseline          = ft_appenddata(cfg,data_hbaseline, data_lbaseline);
        else
            try
                [data_baseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_hb,3,[0.48 3.48]);
            catch
                [data_baseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_lb,7,[0.48 3.48]);
            end
        end

        % Low Bass condition
        if length(unique(data_lowbass.trialinfo))>1
            [data_hlowbass]             = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_hl,4,[0.48 3.48]);
            [data_llowbass]             = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_ll,8,[0.48 3.48]);
            cfg                 = [];
            cfg.keepsampleinfo  ='no';
            data_lowbass          = ft_appenddata(cfg,data_hlowbass, data_llowbass);
        else
            try
                [data_lowbass]              = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_hl,4,[0.48 3.48]);
            catch
                [data_lowbass]              = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_ll,8,[0.48 3.48]);
            end
        end

        % High voice condition
        if length(unique(data_highvoice.trialinfo))>1
            [data_hhighvoice]             = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_hh,5,[0.48 3.48]);
            [data_lhighvoice]             = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_lh,9,[0.48 3.48]);
            cfg                 = [];
            cfg.keepsampleinfo  ='no';
            data_highvoice      = ft_appenddata(cfg,data_hhighvoice, data_lhighvoice);
        else
            try
                [data_highvoice]            = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_hh,5,[0.48 3.48]);
            catch
                [data_highvoice]            = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_lh,9,[0.48 3.48]);

            end
        end

        % Control condition
        if length(unique(data_control.trialinfo))>1
            [data_hcontrol]             = Trinh_MiniTrialMakerAmpEnv(data_control,locs_hc,2,[0.48 3.48]);
            [data_lcontrol]             = Trinh_MiniTrialMakerAmpEnv(data_control,locs_lc,6,[0.48 3.48]);
            cfg                 = [];
            cfg.keepsampleinfo  ='no';
            data_control        = ft_appenddata(cfg,data_hcontrol, data_lcontrol);
        else
            try
                [data_control]              = Trinh_MiniTrialMakerAmpEnv(data_control,locs_hc,2,[0.48 3.48]);
            catch
                [data_control]              = Trinh_MiniTrialMakerAmpEnv(data_control,locs_lc,6,[0.48 3.48]);

            end
        end


        save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\Step2_epochdata\',group{g},'\',baby_name,'_EPOCHS'),...
            'data_baseline','data_control','data_lowbass','data_highvoice');
    end
end


% %% load the amplitude envelopes
% load('D:\MUSICOM_EEG\SideProject-idyom-trf\Stimuli\StimAcoustics_movement.mat')
% 
% hopp_music=denv{1,1};
% hopp_control=denv{3,1};
% hopp_lp=denv{2,1};
% hopp_hp=denv{5,1};
% 
% lola_music=denv{6,1};
% lola_control=denv{8,1};
% lola_lp=denv{3.48,1};
% lola_hp=denv{10,1};
% 
% hopp_music_time = 0:0.04:21.32;
% [pks_hb,locs_hb] = findpeaks(hopp_music,'MinPeakHeight',0.03,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_hb);
% % indices_above_median = pks_hb < median_value;
% % locs_hb = locs_hb(indices_above_median);
% 
% % plot(hopp_music_time,hopp_music,hopp_music_time(locs_hb),pks_hb,'or')
% [pks_lb,locs_lb] = findpeaks(lola_music,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % plot(hopp_music_time,lola_music,hopp_music_time(locs_lb),pks_lb,'or')
% % median_value = median(pks_lb);
% % indices_above_median = pks_lb < median_value;
% % locs_lb = locs_lb(indices_above_median);
% 
% 
% [pks_hc,locs_hc] = findpeaks(hopp_control,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_hc);
% % indices_above_median = pks_hc < median_value;
% % locs_hc = locs_hc(indices_above_median);
% 
% % plot(hopp_music_time,hopp_control,hopp_music_time(locs_hc),pks_hc,'or')
% [pks_lc,locs_lc] = findpeaks(lola_control,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_lc);
% % indices_above_median = pks_lc < median_value;
% % locs_lc = locs_lc(indices_above_median);
% % plot(hopp_music_time,lola_control,hopp_music_time(locs_lc),pks_lc,'or')
% 
% [pks_hh,locs_hh] = findpeaks(hopp_hp,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_hh);
% % indices_above_median = pks_hh < median_value;
% % locs_hh = locs_hh(indices_above_median);
% 
% [pks_lh,locs_lh] = findpeaks(lola_hp,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_lh);
% % indices_above_median = pks_lh < median_value;
% % locs_lh = locs_lh(indices_above_median);
% 
% [pks_hl,locs_hl] = findpeaks(hopp_lp,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_hl);
% % indices_above_median = pks_hl < median_value;
% % locs_hl = locs_hl(indices_above_median);
% 
% [pks_ll,locs_ll] = findpeaks(lola_lp,'MinPeakHeight',0.15,'MinPeakDistance',4,'MinPeakProminence',0.02);
% % median_value = median(pks_ll);
% % indices_above_median = pks_ll < median_value;
% % locs_ll = locs_ll(indices_above_median);
% 
% 
% %% access preprocessed files
% group={'3m','6m','12m'};
% for g=1:length(group)
%     cd(['D:\TinyDancers_Movement\fieldtrip_velo\PM\' group{g}]);
%     FileNames = {dir(['*' '.mat']).name} ;
% 
%     for k=1:1:length(FileNames)
% 
%         % load data
%         load(FileNames{k});
%         baby_name       = erase(FileNames{k},'.mat');
% 
%         % smoothing data
%         data_baseline = Giac_MovingMean( data_baseline, 3 );
%         data_control = Giac_MovingMean( data_control, 3 );
%         data_highvoice = Giac_MovingMean( data_highvoice, 3 );
%         data_lowbass = Giac_MovingMean( data_lowbass, 3 );
% 
%         data_baseline = Giac_MovingMean( data_baseline, 3 );
%         data_control = Giac_MovingMean( data_control, 3 );
%         data_highvoice = Giac_MovingMean( data_highvoice, 3 );
%         data_lowbass = Giac_MovingMean( data_lowbass, 3 );
% 
% 
%         %% Epoching
% 
%         % Reference condition
%         if length(unique(data_baseline.trialinfo))>1
%             [data_hbaseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_hb,3,[0.48 3.48]);
%             [data_lbaseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_lb,3.48,[0.48 3.48]);
%             cfg                 = [];
%             cfg.keepsampleinfo  ='no';
%             data_baseline          = ft_appenddata(cfg,data_hbaseline, data_lbaseline);
%         else
%             try
%                 [data_baseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_hb,3,[0.48 3.48]);
%             catch
%                 [data_baseline]             = Trinh_MiniTrialMakerAmpEnv(data_baseline,locs_lb,3.48,[0.48 3.48]);
%             end
%         end
% 
%         % Low Bass condition
%         if length(unique(data_lowbass.trialinfo))>1
%             [data_hlowbass]             = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_hl,4,[0.48 3.48]);
%             [data_llowbass]             = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_ll,8,[0.48 3.48]);
%             cfg                 = [];
%             cfg.keepsampleinfo  ='no';
%             data_lowbass          = ft_appenddata(cfg,data_hlowbass, data_llowbass);
%         else
%             try
%                 [data_lowbass]              = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_hl,4,[0.48 3.48]);
%             catch
%                 [data_lowbass]              = Trinh_MiniTrialMakerAmpEnv(data_lowbass,locs_ll,8,[0.48 3.48]);
%             end
%         end
% 
%         % High voice condition
%         if length(unique(data_highvoice.trialinfo))>1
%             [data_hhighvoice]             = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_hh,5,[0.48 3.48]);
%             [data_lhighvoice]             = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_lh,9,[0.48 3.48]);
%             cfg                 = [];
%             cfg.keepsampleinfo  ='no';
%             data_highvoice      = ft_appenddata(cfg,data_hhighvoice, data_lhighvoice);
%         else
%             try
%                 [data_highvoice]            = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_hh,5,[0.48 3.48]);
%             catch
%                 [data_highvoice]            = Trinh_MiniTrialMakerAmpEnv(data_highvoice,locs_lh,9,[0.48 3.48]);
% 
%             end
%         end
% 
%         % Control condition
%         if length(unique(data_control.trialinfo))>1
%             [data_hcontrol]             = Trinh_MiniTrialMakerAmpEnv(data_control,locs_hc,2,[0.48 3.48]);
%             [data_lcontrol]             = Trinh_MiniTrialMakerAmpEnv(data_control,locs_lc,6,[0.48 3.48]);
%             cfg                 = [];
%             cfg.keepsampleinfo  ='no';
%             data_control        = ft_appenddata(cfg,data_hcontrol, data_lcontrol);
%         else
%             try
%                 [data_control]              = Trinh_MiniTrialMakerAmpEnv(data_control,locs_hc,2,[0.48 3.48]);
%             catch
%                 [data_control]              = Trinh_MiniTrialMakerAmpEnv(data_control,locs_lc,6,[0.48 3.48]);
% 
%             end
%         end
% 
% 
%         save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\denv\Step2_epochdata\',group{g},'\',baby_name,'_EPOCHS'),...
%             'data_baseline','data_control','data_lowbass','data_highvoice');
%     end
% end
