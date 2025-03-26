function data_out = Trinh_makecontinaveragedData(data_in, conditions,refrainonsets,triallength,prestim)
%% Function for segmenting the data and appending them into continous data stream
% segments the condtitions and appends them to make continous data streams
% insert the trial info for the needed conditions
% Trinh nov22

for i = 1:length(conditions)
    tmp                      = conditions(i);
    
    cond                     = [];
    cond                     = find(data_in.trialinfo(:,1)==tmp);
    
    cfg                      = [];
    cfg.trials               = cond'; 
    data_seg                 = ft_selectdata(cfg,data_in);

    onsets                   = round((refrainonsets* data_seg.fsample) + data_seg.sampleinfo(1) + (prestim* data_seg.fsample));

    cfg                      = [];
    cfg.trl(:,1)             = onsets;% samples
    cfg.trl(:,2)             = onsets+(triallength*data_seg.fsample)-1;% samples
    cfg.trl(:,3)             = 0;% samples
    cfg.trl(:,4)             = conditions(i);
    data_refrain{i}          = ft_redefinetrial(cfg,data_seg);
end

data_multirefrain            = ft_appenddata([],data_refrain{:});

cfg                          = [];
data_out                     = ft_timelockanalysis(cfg, data_multirefrain);
end