%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear all
clc
addpath C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\toolboxes\fieldtrip-20220929
ft_defaults;
addpath(genpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/Giac_ToolBox'))
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/functions')
addpath 'D:/MUSICOM_EEG/onsets'

%% access preprocessed files
% group={'3m','6m','12m'};
group={'adults'};

for g=1:length(group)
    cd(['D:/MUSICOM_EEG/PROC/Step1_processeddatamastoidreref/' group{g}]);
    FileNames = {dir(['*' '_PROC.mat']).name} ;

    for k=1:1:length(FileNames)

        % load data
        load(FileNames{k});
        baby_name       = erase(FileNames{k},'_PROC.mat');
        data_int      = Trinh_filterEEG(data_int,'lp',20);

        %% Epoching
% 
%         % Reference condition
% %         [data_baseline_voice]           = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets','Lola_voice_onsets'},[3,7],[0.5 1.5],3);
        [data_baseline]            = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets','Lola_bass_onsets'},[3,7],[0.5 1.5],3.05);
%         [data_baseline]                 = Trinh_MiniTrialMaker(data_int,{'Hopp_onsets','Lola_onsets'},[3,7],[0.5 1.5],3);
% 
%         % Low Bass condition
% %         [data_lowbass_voice]            = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets','Lola_voice_onsets'},[4,8],[0.5 1.5],3);
        [data_lowbass]             = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets','Lola_bass_onsets'},[4,8],[0.5 1.5],3.05);
%         [data_lowbass]                  = Trinh_MiniTrialMaker(data_int,{'Hopp_onsets','Lola_onsets'},[4,8],[0.5 1.5],3);
% 
%         % High voice condition
% %         [data_highvoice_voice]          = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets','Lola_voice_onsets'},[5,9],[0.5 1.5],3);
        [data_highvoice]           = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets','Lola_bass_onsets'},[5,9],[0.5 1.5],3.05);
%         [data_highvoice]                = Trinh_MiniTrialMaker(data_int,{'Hopp_onsets','Lola_onsets'},[5,9],[0.5 1.5],3);
% 
%         % Control condition
% %         [data_control_voice]            = Trinh_MiniTrialMaker(data_int,{'Hopp_control_voice_onsets','Lola_control_voice_onsets'},[2,6],[0.5 1.5],3);
        [data_control]             = Trinh_MiniTrialMaker(data_int,{'Hopp_control_bass_onsets','Lola_control_bass_onsets'},[2,6],[0.5 1.5],3.05);
%         [data_control]                  = Trinh_MiniTrialMaker(data_int,{'Hopp_control_onsets','Lola_control_onsets'},[2,6],[0.5 1.5],3);
% 
%         % Silence condition
% %         data_silence                = Trinh_recutCondition(data_int,1,2,0.5);
% 
%         %% Clean data for leftover artefacts in each epoch
%         % baseline condition
% %         [data_bas_v_clean,~] = Trinh_automaticartefactdetect(data_baseline_voice, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
% %         [data_bas_b_clean,~] = Trinh_automaticartefactdetect(data_baseline_bass, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
        [data_bas_clean,~] = Trinh_automaticartefactdetect(data_baseline, 'stddev',...
            'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
            'FC4','CP3','CP4','P3','P4'}, 'complete');
%         % control condition
% %         [data_con_v_clean,~] = Trinh_automaticartefactdetect(data_control_voice, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
% %         [data_con_b_clean,~] = Trinh_automaticartefactdetect(data_control_bass, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
        [data_con_clean,~] = Trinh_automaticartefactdetect(data_control, 'stddev',...
            'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
            'FC4','CP3','CP4','P3','P4'}, 'complete');
%         % lowbass condition
% %         [data_lba_v_clean,~] = Trinh_automaticartefactdetect(data_lowbass_voice, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
% %         [data_lba_b_clean,~] = Trinh_automaticartefactdetect(data_lowbass_bass, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
        [data_lba_clean,~] = Trinh_automaticartefactdetect(data_lowbass, 'stddev',...
            'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
            'FC4','CP3','CP4','P3','P4'}, 'complete');
%         % highvoice condition
% %         [data_hvo_v_clean,~] = Trinh_automaticartefactdetect(data_highvoice_voice, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
% %         [data_hvo_b_clean,~] = Trinh_automaticartefactdetect(data_highvoice_bass, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
        [data_hvo_clean,~] = Trinh_automaticartefactdetect(data_highvoice, 'stddev',...
            'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
            'FC4','CP3','CP4','P3','P4'}, 'complete');
%         % silence condition
% %         [data_sil_clean,~] = Trinh_automaticartefactdetect(data_silence, 'stddev',...
% %             'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
% %             'FC4','CP3','CP4','P3','P4'}, 'complete');
% % 
% %         save(strcat('D:\MUSICOM_EEG\PROC\Step2_epochdata_mast\',group{g},'\',baby_name,'_EPOCHS'),...
% %             'data_bas_v_clean','data_con_v_clean','data_lba_v_clean','data_hvo_v_clean','data_sil_clean',...
% %             'data_bas_b_clean','data_con_b_clean','data_lba_b_clean','data_hvo_b_clean');
% 
        save(strcat('D:/MUSICOM_EEG/PROC/Step2_epochdata_mast/',group{g},'/',baby_name,'_EPOCHS'),...
            'data_bas_clean','data_con_clean','data_lba_clean','data_hvo_clean');
% 
        clear data_bas_clean data_con_clean data_lba_clean data_hvo_clean  
%         clear data_bas_v_clean data_con_v_clean data_lba_v_clean data_hvo_v_clean data_sil_clean 
%         clear data_silence data_baseline_voice data_control_voice data_highvoice_voice data_lowbass_voice
%         clear data_bas_b_clean data_con_b_clean data_lba_b_clean data_hvo_b_clean
%         clear data_baseline_bass data_control_bass data_highvoice_bass data_lowbass_bass
        clear data_baseline data_control data_highvoice data_lowbass

                %% prepare data for frequency analysis
        
%                 % Reference condition
%                 data_baseline_contin    = Trinh_makecontinData(data_int,[3,7],0,21,3);
%         
%                 % Control condition
%                 data_control_contin     = Trinh_makecontinData(data_int,[2,6],0,21,3);
%         
%                 % Low bass condition
%                 data_lowbass_contin     = Trinh_makecontinData(data_int,[4,8],0,21,3);
%         
%                 % High voice condition
%                 data_highvoice_contin   = Trinh_makecontinData(data_int,[5,9],0,21,3);
%         
%                 save(strcat('D:\MUSICOM_EEG\PROC\Step2_continuousdata\',group{g},'\',baby_name,'_CONTIN'),...
%                     'data_baseline_contin','data_control_contin','data_lowbass_contin','data_highvoice_contin');
%         
%                 clear data_highvoice_contin data_baseline_contin data_control_contin data_lowbass_contin data_int

        %         savePath = 'D:\MUSICOM_EEG\PROC\Step3_frequencydata_lw_mast';
        %         saveName = [group{g} ' baseline ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( data_baseline_contin, saveName, savePath, [], [] );
        %         saveName = [group{g} ' control ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( data_control_contin, saveName, savePath, [], [] );
        %         saveName = [group{g} ' lowbass ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( data_lowbass_contin, saveName, savePath, [], [] );
        %         saveName = [group{g} ' highvoice ' baby_name];
        %         lwdata = rs_convert_ft2lw_V7( data_highvoice_contin, saveName, savePath, [], [] );


    end
end

