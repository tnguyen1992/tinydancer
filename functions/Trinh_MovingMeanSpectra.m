function [ data_out ] = Trinh_MovingMeanSpectra(data_in, window )
%
% This function computes a moving window of the data (FT format) calling
% another function 'movingmean'. 
%
%% Giacomo Novembre

data_out = data_in;
    for tr = 1:length(data_in.trial)
        data_out.trial{tr} = movingmean( data_in.trial{tr}, window, 2, []); % third element stands for the dimension over which the computation should be performed
    end
end
