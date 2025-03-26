clc
clear all
close all

addpath(['.' filesep 'eeglab2022.0' filesep]);
eeglab;

eeglab_template = EEG;          % empty template of eeglab struct

close all;      % close eeglab

EEG_dyad = pop_biosig('eegdata/dyad02.bdf');     % load one bdf of your dataset

NB_EEG_CHAN = 64;       % set number of EEG electrodes

% Erase '1-' in labels if you did dual EEG
for ch=1:NB_EEG_CHAN
    EEG_dyad.chanlocs(ch).labels = erase(EEG_dyad.chanlocs(ch).labels , '1-');
end

EEG_dyad = pop_chanedit(EEG_dyad);        % then select "BESA" and it'll assign XYZ locs, select OK

EEG_dyad.chanlocs(NB_EEG_CHAN+1:end) = [];       % remove unnecessary electrodes
chanlocs = EEG_dyad.chanlocs;

eeglab_template.chanlocs = chanlocs;

save(['./eeglab_template_biosemi64'],'eeglab_template');     % save it as layout to use in NPA functions