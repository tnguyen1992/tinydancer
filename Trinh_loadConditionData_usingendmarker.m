function data_out = Trinh_loadConditionData_usingendmarker(filename,eventtype,eventvalue,prestim,poststim)
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

% find the event that we're looking for
if length(eventvalue)>1
    for l=1:length(eventvalue)
        startmark=find(strcmp(eventvalue{l},event.value));
        % find all the end markers S 9
        endmark=find(strcmp('S  9',event.value));
        if isempty(endmark)
            endmark=find(strcmp('S137',event.value));
        end
        duration=[];

        for n=1:length(startmark)
            % only take end markers later than event marker
            endmark_relevant=endmark(endmark>startmark(n));
            % take the first end marker after event was started
            % determine duration of that condition
            duration=[duration event.sample(endmark_relevant(1))-event.sample(startmark(n))];
        end

        duration=max(duration)+3;

        cfg                     = [];
        cfg.dataset             = filename;
        cfg.trialfun            = 'ft_trialfun_general';
        cfg.showcallinfo        = 'no';
        cfg.feedback            = 'error';
        cfg.trialdef.eventtype  = eventtype;
        cfg.trialdef.eventvalue = eventvalue{l};
        cfg.trialdef.prestim    = prestim;                                                % get more data so that ASR doesn't mess up
        cfg.trialdef.poststim   = round(duration/1000);                                % convert into seconds

        cfg                     = ft_definetrial(cfg);
        data_out{l}             = ft_redefinetrial(cfg, data_proc);
        
    end
    
    cfg            = [];
    data_all       = ft_appenddata(cfg, data_out{1}, data_out{2});
else
    startmark=find(strcmp(eventvalue,event.value));
    % find all the end markers S 9
    endmark=find(strcmp('S  9',event.value));
    if isempty(endmark)
        endmark=find(strcmp('S137',event.value));
    end

    duration=[];

    for n=1:length(startmark)
        % only take end markers later than event marker
        endmark_relevant=endmark(endmark>startmark(n));
        % take the first end marker after event was started
        % determine duration of that condition
        duration=[duration event.sample(endmark_relevant(1))-event.sample(startmark(n))];
    end

    duration=max(duration)+3;

    cfg                     = [];
    cfg.dataset             = filename;
    cfg.trialfun            = 'ft_trialfun_general';
    cfg.showcallinfo        = 'no';
    cfg.feedback            = 'error';
    cfg.trialdef.eventtype  = eventtype;
    cfg.trialdef.eventvalue = eventvalue;
    cfg.trialdef.prestim    = prestim;                                                % get more data so that ASR doesn't mess up
    cfg.trialdef.poststim   = round(duration/1000);                                % convert into seconds

    cfg                     = ft_definetrial(cfg);
    data_out                = ft_redefinetrial(cfg, data_proc);
end

end