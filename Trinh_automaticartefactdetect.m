function [data_out,trial_count] = Trinh_automaticartefactdetect(data_in, method, sliding, threshold, window, channels,reject_method)
%% function to automatically detect artefacts according to a sliding window approach
% can be used for range or standard deviation 
% method is then either 'range' or 'stddev'
% sliding should be 'yes'
% the default threshold for range is 200 for stddev it's 100
% window length is normally 200ms
% choose relevant channels in cell array
% Trinh Nov 22

%% Automatic artefact detection
cfg = [];
cfg.method                        = method;
cfg.sliding                       = sliding;
cfg.artfctdef.threshold.channel   = channels;                                  % specify channels of interest
cfg.artfctdef.threshold.bpfilter  = 'no';                                   % use no additional bandpass
cfg.artfctdef.threshold.bpfreq    = [];                                     % use no additional bandpass
cfg.artfctdef.threshold.onset     = [];                                     % just defined to get a similar output from ft_artifact_threshold and artifact_threshold
cfg.artfctdef.threshold.offset    = [];                                     % just defined to get a similar output from ft_artifact_threshold and artifact_threshold
cfg.showcallinfo                  = 'no';
cfg.overlap                       = 1;
cfg.artfctdef.threshold.stddev    = threshold;                                 % stddev
cfg.artfctdef.threshold.winsize   = window;
cfg.artfctdef.threshold.trl       = 'all';

cfgAutoArt = artifact_detect(cfg, data_in);

%% plot if wanted
% cfg                                 = [];
% cfg.ylim                            = [-100 100];
% cfg.blocksize                       = 2;
% cfg.viewmode                        = 'vertical';
% cfg.artfctdef.threshold.artifact    = cfgAutoArt_base.artfctdef.threshold.artifact;
% cfg.continuous                      = 'no';
% cfg.channel                         = 'all';
% cfg.allowoverlap                    = 'yes';
% 
% cfgAutoArt_base                     = ft_databrowser(cfg, data_baseline_asr);
% 
% cfg.artfctdef.threshold.artifact    = cfgAutoArt_contr.artfctdef.threshold.artifact;
% cfgAutoArt_contr                    = ft_databrowser(cfg, data_control);
% 
% cfg.artfctdef.threshold.artifact    = cfgAutoArt_lbass.artfctdef.threshold.artifact;
% cfgAutoArt_lbass                    = ft_databrowser(cfg, data_lowbass);
% 
% cfg.artfctdef.threshold.artifact    = cfgAutoArt_hvoice.artfctdef.threshold.artifact;
% cfgAutoArt_hvoice                   = ft_databrowser(cfg, data_highvoice);

%% substitute artefacts with NaNs
cfg                                 = [];
cfg.artfctdef.reject                = reject_method;
cfg.artfctdef.threshold.artifact    = cfgAutoArt.artfctdef.threshold.artifact;
data_out                            = ft_rejectartifact(cfg, data_in);

% % find trials with NaNs
% [ nan_trials_base ]         = Giac_findNanTrials( data_baseline_asr_ch_art, 'OnlyAll' );
% [data_baseline_asr_ch_art_fin]  = Giac_removeTrials(data_baseline_asr_ch_art,nan_trials_base,'reject');

trial_count=size(data_out.trial,2);
