%% Step 3 Tiny Dancers: ERP analyses
clear all
clc

trial_mus=[];
trial_con=[];
trial_hpi=[];
trial_lpi=[];

% load names of all files '3m','6m','12m',
group={'adults'};
for g=1:length(group)
    cd(['D:\MUSICOM_EEG\PROC\Step2_epochdata_mast\' group{g}]);
    FileNames      = {dir(['*' 'EPOCHS.mat']).name} ;

    for k = 1:1:length(FileNames)

        % load data
        load(FileNames{k})
        % get baby ID
        baby_name      = erase(FileNames{k},'_EPOCHS.mat');

        %%
        trial_mus=[trial_mus;g,k,length(data_bas_clean.trial)];
        trial_con=[trial_con;g,k,length(data_con_clean.trial)];
        trial_hpi=[trial_hpi;g,k,length(data_hvo_clean.trial)];
        trial_lpi=[trial_lpi;g,k,length(data_lba_clean.trial)];
    


        %% baseline correction 100 ms
        %         data_bas_v_bc = Trinh_baselinecorrection(data_bas_v_clean,[-0.1 0]);
        %         data_con_v_bc = Trinh_baselinecorrection(data_con_v_clean,[-0.1 0]);
        %         data_hvo_v_bc = Trinh_baselinecorrection(data_hvo_v_clean,[-0.1 0]);
        %         data_lba_v_bc = Trinh_baselinecorrection(data_lba_v_clean,[-0.1 0]);
        %
        %         data_bas_b_bc = Trinh_baselinecorrection(data_bas_b_clean,[-0.1 0]);
        %         data_con_b_bc = Trinh_baselinecorrection(data_con_b_clean,[-0.1 0]);
        %         data_hvo_b_bc = Trinh_baselinecorrection(data_hvo_b_clean,[-0.1 0]);
        %         data_lba_b_bc = Trinh_baselinecorrection(data_lba_b_clean,[-0.1 0]);
        %
        %         data_sil_bc = Trinh_baselinecorrection(data_sil_clean,[-0.1 0]);

        data_bas_bc = Trinh_baselinecorrection(data_bas_clean,[-0.1 0]);
        data_con_bc = Trinh_baselinecorrection(data_con_clean,[-0.1 0]);
        data_hvo_bc = Trinh_baselinecorrection(data_hvo_clean,[-0.1 0]);
        data_lba_bc = Trinh_baselinecorrection(data_lba_clean,[-0.1 0]);
        %% look at ERPs in ft
        cfg = [];
        cfg.keeptrials='no';
%         ERP_baseline_voice   = ft_timelockanalysis(cfg, data_bas_v_bc);
%         ERP_control_voice    = ft_timelockanalysis(cfg, data_con_v_bc);
%         ERP_lowbass_voice    = ft_timelockanalysis(cfg, data_lba_v_bc);
%         ERP_highvoice_voice  = ft_timelockanalysis(cfg, data_hvo_v_bc);
% ERP_baseline_bass    = ft_timelockanalysis(cfg, data_bas_b_bc);
% ERP_control_bass     = ft_timelockanalysis(cfg, data_con_b_bc);
% ERP_lowbass_bass     = ft_timelockanalysis(cfg, data_lba_b_bc);
% ERP_highvoice_bass   = ft_timelockanalysis(cfg, data_hvo_b_bc);
% ERP_music       = ft_timelockanalysis(cfg, data_music);
% ERP_silence          = ft_timelockanalysis(cfg, data_sil_bc);

        ERP_baseline   = ft_timelockanalysis(cfg, data_bas_bc);
        ERP_control    = ft_timelockanalysis(cfg, data_con_bc);
        ERP_lowbass    = ft_timelockanalysis(cfg, data_lba_bc);
        ERP_highvoice  = ft_timelockanalysis(cfg, data_hvo_bc);

%         save(strcat('D:\MUSICOM_EEG\PROC\Step3_ERP_mast\',group{g},'\',baby_name,'_ERP'),...
%             'ERP_baseline_voice','ERP_control_voice','ERP_lowbass_voice','ERP_highvoice_voice',...
%             'ERP_baseline_bass','ERP_control_bass','ERP_lowbass_bass','ERP_highvoice_bass','ERP_silence');

save(strcat('D:\MUSICOM_EEG\PROC\Step3_ERPs\',group{g},'\',baby_name,'_ERP'),...
    'ERP_baseline','ERP_control','ERP_lowbass','ERP_highvoice');

        %     mnepath= strcat('D:\MUSICOM_EEG\PROC\Step3_ERP_mne\',group,'\');
        %     fiff_file = strcat(mnepath,baby_name,'_ERP_baseline-ave.fif');
        %     fieldtrip2fiff(fiff_file, ERP_baseline);
        %     fiff_file = strcat(mnepath,baby_name,'_ERP_control-ave.fif');
        %     fieldtrip2fiff(fiff_file, ERP_control);
        %     fiff_file = strcat(mnepath, baby_name,'_ERP_lowbass-ave.fif');
        %     fieldtrip2fiff(fiff_file, ERP_lowbass);
        %     fiff_file = strcat(mnepath, baby_name,'_ERP_highvoice-ave.fif');
        %     fieldtrip2fiff(fiff_file, ERP_highvoice);
        %     fiff_file = strcat(mnepath, baby_name,'_ERP_music-ave.fif');
        %     fieldtrip2fiff(fiff_file, ERP_music);

        %         cfg = [];
        %         cfg.keeptrials='yes';
        %         ERP_baseline    = ft_timelockanalysis(cfg, data_bas_v_clean);
        %         ERP_control     = ft_timelockanalysis(cfg, data_con_v_clean);
        %         ERP_lowbass     = ft_timelockanalysis(cfg, data_lba_v_clean);
        %         ERP_highvoice   = ft_timelockanalysis(cfg, data_hvo_v_clean);
        % %         ERP_music       = ft_timelockanalysis(cfg, data_music_clean);
        % %         ERP_silence     = ft_timelockanalysis(cfg, data_sil_clean);
        %
        %         savePath = 'D:\MUSICOM_EEG\PROC\Step3_ERP_lw_mast';
        %         saveName = [group{g} ' baseline voice ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_baseline, saveName, savePath, [], [] );
        %         saveName = [group{g} ' control voice ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_control, saveName, savePath, [], [] );
        %         saveName = [group{g} ' lowbass voice ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_lowbass, saveName, savePath, [], [] );
        %         saveName = [group{g} ' highvoice voice ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_highvoice, saveName, savePath, [], [] );
        %         saveName = [group{g} ' music ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_music, saveName, savePath, [], [] );
        %         saveName = [group{g} ' silence ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( ERP_silence, saveName, savePath, [], [] );
    end
end    
       