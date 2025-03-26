function data_out = Trinh_freq_normalize(data_in, mean_freq)
% removes background noise from each frequency bin by subtracting the
% averaged amplitude measured at neighboring frequency bins from each
% frequency bin as seen in Moureaux et al., 2011, Norazadan et al., 2012,
% Cirelli et al., 2019
% Trinh Jan 24
data_out=data_in;
spectrum = data_in.powspctrm;
data_out.powspctrm=(spectrum-mean_freq)./(spectrum+mean_freq);

end

