%% Step 3 Tiny Dancers: ERP analyses
clear all
clc

addpath(genpath('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\toolboxes\Giac_ToolBox'));

% load names of all files
group={'3m','6m','12m'};
for g=1:length(group)
    cd(['D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\Step2_epochdata\' group{g}]);
    FileNames      = {dir(['*' 'EPOCHS.mat']).name} ;

    for k = 1:1:length(FileNames)

        % load data
        load(FileNames{k})
        % get baby ID
        baby_name      = erase(FileNames{k},'_EPOCHS.mat');

        % reject Nan trials before baseline correction
        [ nan_trials ] = Giac_findNanTrials( data_baseline, 'OnlyAll' );
        [data_baseline] = Giac_removeTrials(data_baseline,nan_trials,'reject');

        [ nan_trials ] = Giac_findNanTrials( data_control, 'OnlyAll' );
        [data_control] = Giac_removeTrials(data_control,nan_trials,'reject');

        [ nan_trials ] = Giac_findNanTrials( data_highvoice, 'OnlyAll' );
        [data_highvoice] = Giac_removeTrials(data_highvoice,nan_trials,'reject');

        [ nan_trials ] = Giac_findNanTrials( data_lowbass, 'OnlyAll' );
        [data_lowbass] = Giac_removeTrials(data_lowbass,nan_trials,'reject');


        %% baseline correction 100 ms
%         data_bas_bc = Trinh_baselinecorrection(data_baseline,[-0.08 0]);
%         data_con_bc = Trinh_baselinecorrection(data_control,[-0.08 0]);
%         data_hvo_bc = Trinh_baselinecorrection(data_highvoice,[-0.08 0]);
%         data_lba_bc = Trinh_baselinecorrection(data_lowbass,[-0.08 0]);

        data_bas_bc = data_baseline;
        data_con_bc = data_control;
        data_hvo_bc = data_highvoice;
        data_lba_bc = data_lowbass;
     
        %% rectify
        cfg=[];
        cfg.rectify = 'no';
        data_bas_bc = ft_preprocessing(cfg,data_bas_bc);
        data_con_bc = ft_preprocessing(cfg,data_con_bc);
        data_hvo_bc = ft_preprocessing(cfg,data_hvo_bc);
        data_lba_bc = ft_preprocessing(cfg,data_lba_bc);

        % smoothing data
        data_bas_bc = Giac_MovingMean( data_bas_bc, 3 );
        data_con_bc = Giac_MovingMean( data_con_bc, 3 );
        data_hvo_bc = Giac_MovingMean( data_hvo_bc, 3 );
        data_lba_bc = Giac_MovingMean( data_lba_bc, 3 );

        data_bas_bc = Giac_MovingMean( data_bas_bc, 3 );
        data_con_bc = Giac_MovingMean( data_con_bc, 3 );
        data_hvo_bc = Giac_MovingMean( data_hvo_bc, 3 );
        data_lba_bc = Giac_MovingMean( data_lba_bc, 3 );


        %% apend music condititions
        cfg          = [];
        cfg.keepsampleinfo='no';
        data_music   = ft_appenddata(cfg,data_bas_bc, data_hvo_bc, data_lba_bc) ;
        %% Automatic artefact detection

        cfg = [];
        cfg.showcallinfo                  = 'no';
        cfg.artfctdef.threshold.channel   = {'PM1'};
        cfg.artfctdef.threshold.stddev    = 0.50;                                 % 20%
        cfg.artfctdef.threshold.trl       = 'all';

        autoart_bas = artifact_threshold_rank(cfg, data_bas_bc);
        autoart_con = artifact_threshold_rank(cfg, data_con_bc);
        autoart_hvo = artifact_threshold_rank(cfg, data_hvo_bc);
        autoart_lba = artifact_threshold_rank(cfg, data_lba_bc);

        try
            cfg                                 = [];
            cfg.artfctdef.reject                = 'complete';
            cfg.artfctdef.threshold.artifact    = autoart_bas.artfctdef.threshold.artifact;
            data_bas_clean                      = ft_rejectartifact(cfg, data_bas_bc);

%             cfg = [];
%             cfg.method = 'channel';
%             data_clean = ft_rejectvisual(cfg, data_bas_bc);
%             data_clean = ft_rejectvisual(cfg, data_bas_clean);

            cfg.artfctdef.threshold.artifact    = autoart_con.artfctdef.threshold.artifact;
            data_con_clean                      = ft_rejectartifact(cfg, data_con_bc);

            cfg = [];
            cfg.keeptrials='no';
            cfg.latency            = [-0.1 1];
            ERP_baseline    = ft_timelockanalysis(cfg, data_bas_clean);
            ERP_control     = ft_timelockanalysis(cfg, data_con_clean);

            cfg=[];
            cfg.rectify = 'no';
            ERP_baseline = ft_preprocessing(cfg,ERP_baseline);
            ERP_control = ft_preprocessing(cfg,ERP_control);
       

            numreps=[length(data_bas_clean.trial),length(data_con_clean.trial)];
            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc\',group{g},'\',baby_name,'_ERP'),...
                'ERP_baseline','ERP_control','numreps');%,'ERP_silence');
        end
        try
            cfg                                 = [];
            cfg.artfctdef.reject                = 'complete';
            cfg.artfctdef.threshold.artifact    = autoart_hvo.artfctdef.threshold.artifact;
            data_hvo_clean                            = ft_rejectartifact(cfg, data_hvo_bc);
            cfg.artfctdef.threshold.artifact    = autoart_lba.artfctdef.threshold.artifact;
            data_lba_clean                            = ft_rejectartifact(cfg, data_lba_bc);

            %% look at ERPs in ft
            cfg = [];
            cfg.keeptrials='no';
            cfg.latency            = [-0.1 1];

            ERP_lowbass     = ft_timelockanalysis(cfg, data_hvo_clean);
            ERP_highvoice   = ft_timelockanalysis(cfg, data_lba_clean);

            cfg=[];
            cfg.rectify = 'no';
            ERP_lowbass = ft_preprocessing(cfg,ERP_lowbass);
            ERP_highvoice = ft_preprocessing(cfg,ERP_highvoice);

            numreps=[length(data_hvo_clean.trial),length(data_lba_clean.trial)];

            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_ERPbc_freq\',group{g},'\',baby_name,'_ERP'),...
                'ERP_lowbass','ERP_highvoice', 'numreps');
        end
        try
            % freq mus
            cfg             = [];
            cfg.output      = 'fourier';
            cfg.channel     = {'all'};
            cfg.method      = 'mtmfft';
            cfg.taper       = 'hilbert';
            cfg.foi        = [1.5:0.25:3];
%             cfg.tapsmofrq   = 0.5;

            freq_baseline   = ft_freqanalysis(cfg, data_bas_clean);
            freq_control    = ft_freqanalysis(cfg, data_con_clean);
            [freq_baseline] = ft_freqdescriptives(cfg, freq_baseline);
            [freq_control]  = ft_freqdescriptives(cfg, freq_control);
            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step2b_frequency\',group{g},'\',baby_name,'_freq'),...
                'freq_baseline','freq_control');

            data_bas_autocor = autocorrhythmov(data_bas_clean,25);
            data_con_autocor = autocorrhythmov(data_con_clean,25);
            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor\',group{g},'\',baby_name,'_autocor'),...
                'data_bas_autocor','data_con_autocor');
        end
        try
            % freq
            cfg             = [];
            cfg.output      = 'fourier';
            cfg.channel     = {'all'};
            cfg.method      = 'mtmfft';
            cfg.taper       = 'hilbert';
            cfg.foi         = [1.5:0.25:3];
%             cfg.tapsmofrq   = 0.5;

            freq_highpitch  = ft_freqanalysis(cfg, data_hvo_clean);
            freq_lowpitch   = ft_freqanalysis(cfg, data_lba_clean);
            [freq_lowpitch] = ft_freqdescriptives(cfg, freq_lowpitch);
            [freq_highpitch]= ft_freqdescriptives(cfg, freq_highpitch);

            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step2b_frequency_freq\',group{g},'\',baby_name,'_freq'),...
                'freq_lowpitch','freq_highpitch');

            data_hvo_autocor = autocorrhythmov(data_hvo_clean,25);
            data_lba_autocor = autocorrhythmov(data_lba_clean,25);
            save(strcat('D:\TinyDancers_Movement\fieldtrip_velo\AE_PROC\SD5\Step3_autocor_freq\',group{g},'\',baby_name,'_autocor'),...
                'data_hvo_autocor','data_lba_autocor');
        end

%     for tr=1:length(data_bas_clean.trial)
%         [acf(tr,:), lags(tr,:)] = autocorr(data_lba_clean.trial{1,tr}(1,:));
%     end
%     acf_avg=mean(acf,1);
%     lags_avg=lags(1,:).*0.04;
%     plot(lags_avg,acf_avg)
    end
end
