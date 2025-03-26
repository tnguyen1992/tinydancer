%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear all
clc
group='12m';
cd(['D:\MUSICOM_EEG\PROC\Step1_processeddatamastoidreref\' group]);
%% access preprocessed files
FileNames = {dir(['*' '_PROC.mat']).name} ;

for k=1:1:length(FileNames)

% load data 
load(FileNames{k});
baby_name       = erase(FileNames{k},'_PROC.mat');

%% Start with Hopp
%% Epoching

% Reference condition
[data_baseline_bass]             = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets'},[3],[0.5 1.5],3);
[data_baseline_voice]            = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets'},[3],[0.5 1.5],3);

% Low Bass condition
[data_lowbass_bass]              = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets'},[4],[0.5 1.5],3);
[data_lowbass_voice]             = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets'},[4],[0.5 1.5],3);

% High voice condition
[data_highvoice_bass]            = Trinh_MiniTrialMaker(data_int,{'Hopp_bass_onsets'},[5],[0.5 1.5],3);
[data_highvoice_voice]           = Trinh_MiniTrialMaker(data_int,{'Hopp_voice_onsets'},[5],[0.5 1.5],3);

% Control condition
[data_control_bass]              = Trinh_MiniTrialMaker(data_int,{'Hopp_control_bass_onsets'},[2],[0.5 1.5],3); 
[data_control_voice]             = Trinh_MiniTrialMaker(data_int,{'Hopp_control_voice_onsets'},[2],[0.5 1.5],3); 

%% Clean data for leftover artefacts in each epoch
% baseline condition
[data_bas_bass_clean,~] = Trinh_automaticartefactdetect(data_baseline_bass, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');
[data_bas_voice_clean,~] = Trinh_automaticartefactdetect(data_baseline_voice, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');

% control condition
[data_con_bass_clean,~] = Trinh_automaticartefactdetect(data_control_bass, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');
[data_con_voice_clean,~] = Trinh_automaticartefactdetect(data_control_voice, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');

% lowbass condition
[data_lba_bass_clean,~] = Trinh_automaticartefactdetect(data_lowbass_bass, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');
[data_lba_voice_clean,~] = Trinh_automaticartefactdetect(data_lowbass_voice, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');

% highvoice condition
[data_hvo_bass_clean,~] = Trinh_automaticartefactdetect(data_highvoice_bass, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');
[data_hvo_voice_clean,~] = Trinh_automaticartefactdetect(data_highvoice_voice, 'stddev',...
    'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
    'FC4','CP3','CP4','P3','P4'}, 'complete');

save(strcat('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\Data\PROC\Step2_epochdatabysongbyvoice_reref\Hopp\',group,'\',baby_name,'_EPOCHS'),...
    'data_bas_bass_clean','data_con_bass_clean','data_lba_bass_clean','data_hvo_bass_clean',...
    'data_bas_voice_clean','data_con_voice_clean','data_lba_voice_clean','data_hvo_voice_clean');

%% prepare data for frequency analysis

% Reference condition
data_baseline_contin    = Trinh_makecontinData(data_int,[3,7],[-1,6,13],8,3);

% Control condition
data_control_contin     = Trinh_makecontinData(data_int,[2,6],[-1,6,13],8,3);

% Low bass condition
data_lowbass_contin     = Trinh_makecontinData(data_int,[4,8],[-1,6,13],8,3);

% High voice condition
data_highvoice_contin   = Trinh_makecontinData(data_int,[5,9],[-1,6,13],8,3);

save(strcat('D:\MUSICOM_EEG\PROC\Step2_continuousdata_avg\',group,'\',baby_name,'_CONTIN'),...
    'data_baseline_contin','data_control_contin','data_lowbass_contin','data_highvoice_contin');


%% then we do the same for lola
%% Epoching
% 
% % Reference condition
% [data_baseline_bass]             = Trinh_MiniTrialMaker(data_int,{'Lola_bass_onsets'},[3],[0.5 1.5],3);
% [data_baseline_voice]            = Trinh_MiniTrialMaker(data_int,{'Lola_voice_onsets'},[3],[0.5 1.5],3);
% 
% % Low Bass condition
% [data_lowbass_bass]              = Trinh_MiniTrialMaker(data_int,{'Lola_bass_onsets'},[4],[0.5 1.5],3);
% [data_lowbass_voice]             = Trinh_MiniTrialMaker(data_int,{'Lola_voice_onsets'},[4],[0.5 1.5],3);
% 
% % High voice condition
% [data_highvoice_bass]            = Trinh_MiniTrialMaker(data_int,{'Lola_bass_onsets'},[5],[0.5 1.5],3);
% [data_highvoice_voice]           = Trinh_MiniTrialMaker(data_int,{'Lola_voice_onsets'},[5],[0.5 1.5],3);
% 
% % Control condition
% [data_control_bass]              = Trinh_MiniTrialMaker(data_int,{'Lola_control_bass_onsets'},[2],[0.5 1.5],3); 
% [data_control_voice]             = Trinh_MiniTrialMaker(data_int,{'Lola_control_voice_onsets'},[2],[0.5 1.5],3); 
% 
% %% Clean data for leftover artefacts in each epoch
% % baseline condition
% [data_bas_bass_clean,~] = Trinh_automaticartefactdetect(data_baseline_bass, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% [data_bas_voice_clean,~] = Trinh_automaticartefactdetect(data_baseline_voice, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% 
% % control condition
% [data_con_bass_clean,~] = Trinh_automaticartefactdetect(data_control_bass, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% [data_con_voice_clean,~] = Trinh_automaticartefactdetect(data_control_voice, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% 
% % lowbass condition
% [data_lba_bass_clean,~] = Trinh_automaticartefactdetect(data_lowbass_bass, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% [data_lba_voice_clean,~] = Trinh_automaticartefactdetect(data_lowbass_voice, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% 
% % highvoice condition
% [data_hvo_bass_clean,~] = Trinh_automaticartefactdetect(data_highvoice_bass, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% [data_hvo_voice_clean,~] = Trinh_automaticartefactdetect(data_highvoice_voice, 'stddev',...
%     'yes', 100, 200, {'Cz','Fz','FCz','Pz','C3','C4','F3','F4','FC3',...
%     'FC4','CP3','CP4','P3','P4'}, 'complete');
% 
% save(strcat('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\Data\PROC\Step2_epochdatabysongbyvoice_reref\Lola\',group,'\',baby_name,'_EPOCHS'),...
%     'data_bas_bass_clean','data_con_bass_clean','data_lba_bass_clean','data_hvo_bass_clean',...
%     'data_bas_voice_clean','data_con_voice_clean','data_lba_voice_clean','data_hvo_voice_clean');
% 
% %% prepare data for frequency analysis
% 
% % Reference condition
% data_baseline_contin    = Trinh_makecontinData(data_int, 7,[-1,6,13],8,3);
% 
% % Control condition
% data_control_contin     = Trinh_makecontinData(data_int, 6,[-1,6,13],8,3);
% 
% % Low bass condition
% data_lowbass_contin     = Trinh_makecontinData(data_int, 8,[-1,6,13],8,3);
% 
% % High voice condition
% data_highvoice_contin   = Trinh_makecontinData(data_int,9,[-1,6,13],8,3);
% 
% save(strcat('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\Data\PROC\Step2_continousdatabysong_reref\Lola\',group,'\',baby_name,'_CONTIN'),...
%     'data_baseline_contin','data_control_contin','data_lowbass_contin','data_highvoice_contin');

end


