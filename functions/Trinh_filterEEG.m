function data_out = Trinh_filterEEG(data_in,filter,freq)
% filters data according to need, can be bandpass, noise, low or high
% insert data from ft_processing
% insert filter (bp,lp,hp) as string
% insert freq in array [0.3 30]
% trinh nov 22

if filter=='bp'
    cfg                     = [];
    cfg.bpfilter            = 'yes';
    cfg.bpfreq              = freq;
    if freq(2)>50
        cfg.dftfilter       = 'yes';
    end
    cfg.bpfilttype          = 'fir';
    data_out                = ft_preprocessing(cfg, data_in);
elseif filter=='lp'
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = freq;
    cfg.lpfiltord           = 5;
    data_out                = ft_preprocessing(cfg, data_in);
elseif filter=='hp'
    cfg                     = [];
    cfg.hpfilter            = 'yes';
    cfg.hpfreq              = freq;
    cfg.hpfiltord           = 5;
    data_out                = ft_preprocessing(cfg, data_in);

end
end
