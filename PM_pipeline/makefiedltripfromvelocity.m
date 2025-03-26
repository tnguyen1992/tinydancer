% go to data folder
% cd(['D:\TinyDancers_Movement\DLC_output\velo_timeseries\']);
cd(['D:\TinyDancers_Movement\DLC_output\PM_postimeseries\']);

addpath('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\toolboxes\fieldtrip-20220929');
ft_defaults;

% load names of all files
% get all babies
FileNames_all      = {dir(['*' '.csv']).name} ;
Names = cellfun(@(x) x(16:19), FileNames_all, 'UniformOutput', false);
% Names = cellfun(@(x) x(1:4), FileNames_all, 'UniformOutput', false);
Names = unique(Names);

% define all conditions
Conditions = {'Silence','HControl','HBaseline','HLowBass', 'HHighVoice',...
    'LControl','LBaseline','LLowBass', 'LHighVoice'};


% find the groups for all the babies
% Define the three folder names
folderNames = {'D:\TinyDancers_Movement\fieldtrip_velo\3m', 'D:\TinyDancers_Movement\fieldtrip_velo\6m', 'D:\TinyDancers_Movement\fieldtrip_velo\12m'};

% Initialize a cell array to store the file names and their folder names
fileList = {};

% Loop through each folder name
for i = 1:length(folderNames)

    % Get a list of all the files in the current folder
    files = dir(folderNames{i});

    % Loop through each file and add its name and folder name to the fileList cell array
    for j = 1:length(files)
        % Remove the file extension from the filename
        [pathstr, name, ext] = fileparts(files(j).name);
        fileList{end+1,1} = name;
        fileList{end,2} = folderNames{i};
    end

end

%%

for j=1:length(Names)
    data_hcontrol=[];
    data_hbaseline=[];
    data_hlowbass=[];
    data_hhighvoice=[];
    data_lcontrol=[];
    data_lbaseline=[];
    data_llowbass=[];
    data_lhighvoice=[];

    try
        data_hcontrol = velo2ft(Names{j},Conditions{2});
    end
    try
        data_hbaseline = velo2ft(Names{j},Conditions{3});
    end
    try
        data_hlowbass = velo2ft(Names{j},Conditions{4});
    end
    try
        data_hhighvoice = velo2ft(Names{j},Conditions{5});
    end
    try
        data_lcontrol = velo2ft(Names{j},Conditions{6});
    end
    try
        data_lbaseline = velo2ft(Names{j},Conditions{7});
    end
    try
        data_llowbass = velo2ft(Names{j},Conditions{8});
    end
    try
        data_lhighvoice = velo2ft(Names{j},Conditions{9});
    end

    % create a structure for each condition
    cfg=[];
    cfg.padtype  = 'no';
    cfg.parameter = {'trial'};
    
    if ~isempty(data_hcontrol) && ~isempty(data_lcontrol)
        % Both data_hcontrol and data_lcontrol exist
        data_control = ft_appenddata(cfg, data_hcontrol, data_lcontrol);
        data_control.trialinfo = [data_hcontrol.trialinfo; data_lcontrol.trialinfo];
    elseif ~isempty(data_hcontrol)
        % Only data_hcontrol exists
        data_control = data_hcontrol;
        data_control.trialinfo = data_hcontrol.trialinfo;
    elseif ~isempty(data_lcontrol)
        % Only data_lcontrol exists
        data_control = data_lcontrol;
        data_control.trialinfo = data_lcontrol.trialinfo;
    end
    
    if ~isempty(data_hbaseline) && ~isempty(data_lbaseline)
        data_baseline = ft_appenddata(cfg,data_hbaseline, data_lbaseline);
        data_baseline.trialinfo = [data_hbaseline.trialinfo;data_lbaseline.trialinfo];
    elseif ~isempty(data_hbaseline)
        data_baseline = data_hbaseline;
        data_baseline.trialinfo = data_hbaseline.trialinfo;
    elseif ~isempty(data_lbaseline)
        % Only data_lcontrol exists
        data_baseline = data_lbaseline;
        data_baseline.trialinfo = data_lbaseline.trialinfo;
    end

    if ~isempty(data_hlowbass) && ~isempty(data_llowbass) 
    data_lowbass = ft_appenddata(cfg,data_hlowbass,data_llowbass);
    data_lowbass.trialinfo = [data_hlowbass.trialinfo;data_llowbass.trialinfo];
    elseif ~isempty(data_hlowbass)
        data_lowbass = data_hlowbass;
        data_lowbass.trialinfo = data_hlowbass.trialinfo;
    elseif ~isempty(data_llowbass)
        % Only data_lcontrol exists
        data_lowbass = data_llowbass;
        data_lowbass.trialinfo = data_llowbass.trialinfo;
    end

    if ~isempty(data_hhighvoice) && ~isempty(data_lhighvoice) 
    data_highvoice = ft_appenddata(cfg,data_hhighvoice,data_lhighvoice);
    data_highvoice.trialinfo = [data_hhighvoice.trialinfo;data_lhighvoice.trialinfo];
    elseif ~isempty(data_hhighvoice)
        data_highvoice = data_hhighvoice;
        data_highvoice.trialinfo = data_hhighvoice.trialinfo;
    elseif ~isempty(data_lhighvoice)
        % Only data_lcontrol exists
        data_highvoice = data_lhighvoice;
        data_highvoice.trialinfo = data_lhighvoice.trialinfo;
    end

    % Get the folder name for the current result based on its name
    idx = find(strcmp(fileList(:,1), Names{j}), 1);
    if ~isempty(idx)
        folderName = fileList{idx, 2};
    else
        error('Result name not found in file list.');
    end

    save([folderName '\' Names{j} '.mat'],...
        "data_baseline","data_control","data_lowbass","data_highvoice")

    clear data_lcontrol data_lbaseline data_llowbass data_lhighvoice data_hhighvoice ...
        data_hcontrol data_hbaseline data_hlowbass
end



