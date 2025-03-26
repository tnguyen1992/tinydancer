
function ft_struct = Trinh_happe_detectbadchannels(ft_struct,baby,flatline,corrnoise,specrej)
% rejects channels according to the HAPPILEE pipeline
% baby (1 on, 0 off) --> changes channel locations according to MUSICOM WP4
% (needs adaption for Biosemi as the function needs 3d coordinates). This
% can be achieved by saving the chanlocs file from eeglab (pop_chanedit)
% flatline (1 on, 0 off) --> default 5 seconds
% corrnoise (1 on, 0 off) --> determine correlation degree (0.7 default),
% determine noise_threshold (default 2.5)
% specrej (1 on, 0 off) --> rejects based on outliers in power in the signal over all
% frequency bins (default specthreshold [-2.75 2.75])
%
% example: % ft_struct =
% Trinh_happe_detectbadchannels(ft_struct,1,1,[],1,[],[],1, []) for default
% mode for WP4 data
% 24/10/22 Trinh (IIT)

%%%% CONCATENATE ALL YOUR TRIALS %%%%
all_trials = cat(2,ft_struct.trial{:});

%%% load template to create eeglab structure %%%
% Update channel locations for MUSICOM 4
if baby
load('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\scripts\Trinh_MUSICOM_eeglab_chanlocs.mat');
eeglab_struct=ans;
eeglab_struct.chanlocs=ans.chanlocs;
else
eeglab_template = load('template_ft_eeglab.mat','template_ft_eeglab');
eeglab_struct = eeglab_template.template_ft_eeglab;
end
%%% update the data field (Nchan x Ntime) %%%
eeglab_struct.data = all_trials;

%%% update other fields of interest %%%
[eeglab_struct.nbchan , eeglab_struct.pnts] = size(all_trials);         % number of channels , number of time points
eeglab_struct.trials = 1;           % because 1 big trial (of concatenated data)
eeglab_struct.srate = ft_struct.fsample;         % important if you have a sampling rate ≠ 2048Hz

% The following part is maybe not necessary (it updates the time
% information, saying "tStart = tStart of the original data in trial 1" and
% then deducting tStop based on srate and the num of time points concatenated
% (eg, you concatenated 3 trials [-1s, 10s], it's gonna create [-1s,32s]))
eeglab_struct.xmin = ft_struct.time{1}(1);
eeglab_struct.xmax =  eeglab_struct.xmin + (eeglab_struct.pnts-1)/eeglab_struct.srate;
eeglab_struct.times = eeglab_struct.xmin : 1/eeglab_struct.srate : eeglab_struct.xmax;



if flatline 
    eeglab_struct = pop_clean_rawdata(eeglab_struct, 'FlatlineCriterion', 5, ...
        'ChannelCriterion', .1, 'LineNoiseCriterion', ...
        20, 'Highpass', 'off', 'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', 'BurstRejection', 'off', ...
        'Distance', 'Euclidian') ;
%     eeglab_struct                    = clean_flatlines(eeglab_struct,max_flatline_duration);
end

if corrnoise
     eeglab_struct = pop_clean_rawdata(eeglab_struct, 'FlatlineCriterion', ...
             'off', 'ChannelCriterion', .7, 'LineNoiseCriterion', ...
             2.5, 'Highpass', 'off', 'BurstCriterion', 'off', ...
             'WindowCriterion', 'off', 'BurstRejection', ...
             'off', 'Distance', 'Euclidian') ;
end

if specrej

    eeglab_struct = pop_rejchan(eeglab_struct, 'elec', [1:eeglab_struct.nbchan], ...
            'threshold', [-2.75 2.75], 'norm', 'on', 'measure', ...
            'spec', 'freqrange', [0.3 30]) ;
end 

ft_struct.trial = mat2cell(eeglab_struct.data,eeglab_struct.nbchan,cell2mat(cellfun(@(x) size(x,2),ft_struct.trial,'UniformOutput',0)));
ft_struct.label = {eeglab_struct.chanlocs.labels};
% reformat EEG data structure back to ft (uneven trials allowed)
end
