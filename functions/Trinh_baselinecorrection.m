function data_out = Trinh_baselinecorrection(data_in,baseline)

cfg = [];
cfg.demean             = 'yes';
cfg.baselinewindow     = baseline; % [begin end] (default = 'no')
cfg.channel            = 'all'; % cell-array, see FT_CHANNELSELECTION
cfg.parameter          = 'trial'; %field for whic

[data_out]             = ft_preprocessing(cfg, data_in);
end
