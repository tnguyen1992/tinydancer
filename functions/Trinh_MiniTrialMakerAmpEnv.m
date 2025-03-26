function [all_erp_all_midi] = Trinh_MiniTrialMakerAmpEnv(data,note_onsets,code_corresponding_to_midi,n_offsets)
%
% midi_names = cell with strings containing the names of the files of
% interest
% code_corresponding_to_midi = trigger codes, present in the FieldTrip
% data structure, associated to the midi_names (Note: the order has to be
% consistent!)
% dir_midi_files = directory where to find the midi files
% prestum = duration of prestimulus time in s
% n_offsets = 2 values, how much time you want to have before and after the
% note of interest (in s)
%% Trinh (based on Giac & Robi's Mini Trial Maker function) Nov 22

addpath('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\IIT_Postdoc\WP4\Stimuli\midi\matlab-midi-master\src\'); % tools for getting midi information

tr_info = data.trialinfo;
n_offsets = round(n_offsets * data.fsample);

all_tr_all_midi = [];

for i = 1:length(code_corresponding_to_midi)

    tmp = code_corresponding_to_midi(i);
    instances_of_same_midi = find(tr_info(:,1)==tmp);

    cfg              = [];
    cfg.trials       = instances_of_same_midi';
    EEG_of_same_midi = ft_selectdata(cfg,data);

    % Trinh: access onsets from saved files (because of combination of two
    % voices)
%     load(midi_names{i});
    all_tr_same_midi = [];

    for ii = 1:length(instances_of_same_midi) % loop into same exposure of 1 song

        cfg                      = [];
        cfg.trials               = ii;
        EEG_of_one_midi          = ft_selectdata(cfg,EEG_of_same_midi);
        % correct note onsets by sampleinfo and prestim
%         try
            note_onsets_corrected = note_onsets  + EEG_of_one_midi.sampleinfo(1) ;
%         catch
%             note_onsets_corrected = note_onsets;
%         end
        cfg                              = [];
        cfg.trl(:,1)                     = note_onsets_corrected-n_offsets(1);% samples
        cfg.trl(:,2)                     = note_onsets_corrected+n_offsets(2)-1;% samples
        cfg.trl(:,3)                     = -n_offsets(1);% samples
        cfg.trl(:,4)                     = tmp;
        EEG_epochs_of_one_midi           = ft_redefinetrial(cfg,EEG_of_one_midi);

        all_tr_same_midi{ii} = EEG_epochs_of_one_midi;

    end

    all_erp_same_midi = ft_appenddata([],all_tr_same_midi{:});
    all_erp_same_midi.fsample = data.fsample;

    all_tr_all_midi{i} = all_erp_same_midi;

end % loop

all_erp_all_midi = ft_appenddata([],all_tr_all_midi{:});

if all_erp_all_midi.fsample==1024
    cfg             = [];
    cfg.resamplefs  = 1000;
    cfg.sampleindex = 'yes';
    [all_erp_all_midi]       = ft_resampledata(cfg, all_erp_all_midi);
end

end % function
