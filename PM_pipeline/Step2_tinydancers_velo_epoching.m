%% Step 2 Tiny Dancers: Epoching and artefact rejection
clear all
clc

addpath 'C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\functions'
addpath 'D:\MUSICOM_EEG/onsets'
addpath('D:\MUSICOM_EEG/EEG/')
addpath('C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/Documenti/MATLAB/toolboxes/fieldtrip-20220929/')
ft_defaults;
%% access preprocessed files
group={'3m','6m','12m'};
for g=1:length(group)
    cd(['D:\TinyDancers_Movement\fieldtrip_velo\PM\' group{g}]);
    FileNames = {dir(['*' '.mat']).name} ;

    for k=1:1:length(FileNames)

        % load data
        load(FileNames{k});
        baby_name       = erase(FileNames{k},'.mat');

%         cfg=[];
%         cfg.rectify = 'yes';
%         data_baseline = ft_preprocessing(cfg,data_baseline);
%         data_control = ft_preprocessing(cfg,data_control);
%         data_highvoice = ft_preprocessing(cfg,data_highvoice);
%         data_lowbass = ft_preprocessing(cfg,data_lowbass);
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
%         data_baseline      = Trinh_filterEEG(data_baseline,'lp',6);
%         data_control      = Trinh_filterEEG(data_control,'lp',6);
%         data_highvoice      = Trinh_filterEEG(data_highvoice,'lp',6);
%         data_lowbass      = Trinh_filterEEG(data_lowbass,'lp',6);

        %% Epoching

        % Reference condition
        if length(unique(data_baseline.trialinfo))>1
            [data_baseline]             = Trinh_MiniTrialMaker(data_baseline,{'Hopp_onsets','Lola_onsets'},[3,7],[0.5 1.5],[]);
        else
            try
                [data_baseline]             = Trinh_MiniTrialMaker(data_baseline,{'Hopp_onsets'},3,[0.5 1.5],[]);
            catch
                [data_baseline]             = Trinh_MiniTrialMaker(data_baseline,{'Lola_onsets'},[7],[0.5 1.5],[]);
            end
        end

        % Low Bass condition
        if length(unique(data_lowbass.trialinfo))>1
            [data_lowbass]              = Trinh_MiniTrialMaker(data_lowbass,{'Hopp_onsets','Lola_onsets'},[4,8],[0.5 1.5],[]);
        else
            try
                [data_lowbass]              = Trinh_MiniTrialMaker(data_lowbass,{'Hopp_onsets'},4,[0.5 1.5],[]);
            catch
                [data_lowbass]              = Trinh_MiniTrialMaker(data_lowbass,{'Lola_onsets'},8,[0.5 1.5],[]);
            end
        end

        % High voice condition
        if length(unique(data_highvoice.trialinfo))>1
            [data_highvoice]            = Trinh_MiniTrialMaker(data_highvoice,{'Hopp_onsets','Lola_onsets'},[5,9],[0.5 1.5],[]);
        else
            try
                [data_highvoice]            = Trinh_MiniTrialMaker(data_highvoice,{'Hopp_onsets'},5,[0.5 1.5],[]);
            catch
                [data_highvoice]            = Trinh_MiniTrialMaker(data_highvoice,{'Lola_onsets'},[9],[0.5 1.5],[]);

            end
        end

        % Control condition
        if length(unique(data_control.trialinfo))>1
            [data_control]              = Trinh_MiniTrialMaker(data_control,{'Hopp_control_onsets','Lola_control_onsets'},[2,6],[0.5 1.5],[]);
        else
            try
                [data_control]              = Trinh_MiniTrialMaker(data_control,{'Hopp_control_onsets'},2,[0.5 1.5],[]);
            catch
                [data_control]              = Trinh_MiniTrialMaker(data_control,{'Lola_control_onsets'},6,[0.5 1.5],[]);

            end
        end


        save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\PM_PROC\Step2_epochdata\',group{g},'\',baby_name,'_EPOCHS'),...
            'data_baseline','data_control','data_lowbass','data_highvoice');
    end
end
% %% prepare data for frequency analysis
%
% Reference condition
% data_baseline_contin    = Trinh_makecontinData(data_int,[3,7],0,21,3);
% 
% % Control condition
% data_control_contin     = Trinh_makecontinData(data_int,[2,6],0,21,3);
% 
% % Low bass condition
% data_lowcontin     = Trinh_makecontinData(data_int,[4,8],0,21,3);
% 
% % High voice condition
% data_highvoice_contin   = Trinh_makecontinData(data_int,[5,9],0,21,3);
% 
% save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\PROC\Step2_continuousdata\',group,'\',baby_name,'_CONTIN'),...
%     'data_baseline_contin','data_control_contin','data_lowcontin','data_highvoice_contin');
%
% end
%

