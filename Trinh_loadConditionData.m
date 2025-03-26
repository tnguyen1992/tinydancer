function data_out = Trinh_loadConditionData(filename,eventtype,eventvalue,prestim,poststim)
%% load and cut condition data
% insert name of the file: filename, e.g., 'MUSICOM_001.vhdr'
% eventtype, e.g., 'Stimulus' or 'Trigger'
% eventvalue, e.g., 'S  1' in cell
% define how much to cut before and after trigger
% Trinh (nov 22)

% load the data
cfg                     = [];
cfg.dataset             = filename;
data_proc               = ft_preprocessing(cfg);

% extract event structure
event = ft_read_event(filename);
event = struct2table(event);

% startmark=find(strcmp(eventvalue,event.value));
% endmark=find(strcmp('S  9',event.value));
% endmark_relevant=endmark(endmark>startmark(1));
% endmark_relevant=endmark_relevant(1);
cfg                     = [];
cfg.dataset             = filename;
cfg.trialfun            = 'ft_trialfun_general';
cfg.showcallinfo        = 'no';
cfg.feedback            = 'error';
cfg.trialdef.eventtype  = eventtype;
cfg.trialdef.eventvalue = eventvalue;
cfg.trialdef.prestim    = prestim;                                                % get more data so that ASR doesn't mess up 
cfg.trialdef.poststim   = poststim;
cfg.representation      ='table';

cfg                     = ft_definetrial(cfg);  
data_out                = ft_redefinetrial(cfg, data_proc);

end