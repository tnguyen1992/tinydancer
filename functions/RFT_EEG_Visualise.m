function [ cfgVISUAL ]= RFT_EEG_Visualise( data, channels, range, layout_file)
%


% PREPARE LAYOUT
load(layout_file); % this normally load a 'lay' structure
cfg.layout          = lay;
layout              = ft_prepare_layout(cfg, data);


% VISUALIZE DATA
cfgint                 = [];
cfgint.blocksize = 10;
cfgint.layout          = layout;
cfgint.ylim            = range;
cfgint.continuous      = 'no';
cfgint.selectmode      = 'markartifact';
cfgint.channel         = channels;
cfgint.viewmode        = 'vertical';
cfgint.axisfontsize    = 10; 
cfgint.plotlabels      = 'yes';
cfgint.artifact        = [];

cfgVISUAL =[];
%cfgVISUAL              = ft_databrowser(cfgint,data);
ft_databrowser(cfgint,data);
