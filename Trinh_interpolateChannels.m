function data_out  = Trinh_interpolateChannels(data_in,missingchannels,neighboursfile,method)     

if ~isempty(missingchannels)
    cfg                = [];
    cfg.method         = method;
    cfg.missingchannel = missingchannels;  
    cfg.senstype       = 'eeg';
    cfg.neighbours     = neighboursfile;
    cfg.trials         = 'all';
    cfg.elec           = ft_read_sens('standard_1020.elc'); % Trinh trick
    [data_out] = ft_channelrepair(cfg, data_in);
else
    data_out=data_in;
end
end