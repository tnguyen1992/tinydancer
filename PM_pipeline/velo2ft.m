function newdata=velo2ft(filename,condition)
% this function imports the velocity data from the directory
% and extracts the right condition
% defines the sampling rate of the data
% and the channels
% define body parts
data_template.label{1, 1}  = 'PM1';
data_template.label{2, 1}  = 'PM2';
data_template.label{3, 1}  = 'PM3';
data_template.label{4, 1}  = 'PM4';
data_template.label{5, 1}  = 'PM5';
data_template.label{6, 1}  = 'PM6';
data_template.label{7, 1}  = 'PM7';
data_template.label{8, 1}  = 'PM8';
data_template.label{9, 1}  = 'PM9';
data_template.label{10, 1}  = 'PM10';
% data_template.label{1, 1}  = 'righteye_x';
% data_template.label{2, 1}  = 'righteye_y';
% data_template.label{3, 1}  = 'lefteye_x';
% data_template.label{4, 1}  = 'lefteye_y';
% data_template.label{5, 1}  = 'mouth_x';
% data_template.label{6, 1}  = 'mouth_y';
% data_template.label{7, 1}  = 'rightshoulder_x';
% data_template.label{8, 1}  = 'rightshoulder_y';
% data_template.label{9, 1}  = 'leftshoulder_x';
% data_template.label{10, 1}  = 'leftshoulder_y';
% data_template.label{11, 1}  = 'chest_x';
% data_template.label{12, 1}  = 'chest_y';
% data_template.label{13, 1}  = 'rightelbow_x';
% data_template.label{14, 1}  = 'rightelbow_y';
% data_template.label{15, 1}  = 'rightwrist_x';
% data_template.label{16, 1}  = 'rightwrist_y';
% data_template.label{17, 1}  = 'righthand_x';
% data_template.label{18, 1}  = 'righthand_y';
% data_template.label{19, 1}  = 'leftelbow_x';
% data_template.label{20, 1}  = 'leftelbow_y';
% data_template.label{21, 1}  = 'leftwrist_x';
% data_template.label{22, 1}  = 'leftwrist_y';
% data_template.label{23, 1}  = 'lefthand_x';
% data_template.label{24, 1}  = 'lefthand_y';
% data_template.label{25, 1}  = 'rightknee_x';
% data_template.label{26, 1}  = 'rightknee_y';
% data_template.label{27, 1}  = 'rightankle_x';
% data_template.label{28, 1}  = 'rightankle_y';
% data_template.label{29, 1}  = 'rightfoot_x';
% data_template.label{30, 1}  = 'rightfoot_y';
% data_template.label{31, 1}  = 'leftknee_x';
% data_template.label{32, 1}  = 'leftknee_y';
% data_template.label{33, 1}  = 'leftankle_x';
% data_template.label{34, 1}  = 'leftankle_y';
% data_template.label{35, 1}  = 'leftfoot_x';
% data_template.label{36, 1}  = 'leftfoot_y';

% data_template.label{1, 1}  = 'righteye';
% data_template.label{2, 1}  = 'lefteye';
% data_template.label{3, 1}  = 'mouth';
% data_template.label{4, 1}  = 'rightshoulder';
% data_template.label{5, 1}  = 'leftshoulder';
% data_template.label{6, 1}  = 'chest';
% data_template.label{7, 1}  = 'rightelbow';
% data_template.label{8, 1}  = 'rightwrist';
% data_template.label{9, 1}  = 'righthand';
% data_template.label{10, 1}  = 'leftelbow';
% data_template.label{11, 1}  = 'leftwrist';
% data_template.label{12, 1}  = 'lefthand';
% data_template.label{13, 1}  = 'rightknee';
% data_template.label{14, 1}  = 'rightankle';
% data_template.label{15, 1}  = 'rightfoot';
% data_template.label{16, 1}  = 'leftknee';
% data_template.label{17, 1}  = 'leftankle';
% data_template.label{18, 1}  = 'leftfoot';

filename_pattern = ['pm-timeseries__' filename '_' condition '*_coordinates.csv.csv'];
% filename_pattern = [filename '_' condition '*_velocity.csv'];

matching_files = dir(fullfile('D:\TinyDancers_Movement\DLC_output\PM_timeseries',filename_pattern));

% Extract only the names of the matching files
file_names = {matching_files.name};

for i=1:length(file_names)
    data.fsample = 25;
    data.label=data_template.label;
    data_in = importpm(file_names{i});
    velo_data = table2array(data_in)';
    data.trial{1, i}=velo_data;

    data.time{1, i} = 0.04:1/data.fsample:length(data.trial{1, i})/data.fsample;
    if strcmp(condition,'Silence')
        data.trialinfo(i,1) = 1;
    elseif strcmp(condition,'HControl')
        data.trialinfo(i,1) = 2;
    elseif strcmp(condition,'HBaseline')
        data.trialinfo(i,1) = 3;
    elseif strcmp(condition,'HLowBass')
        data.trialinfo(i,1) = 4;
    elseif strcmp(condition,'HHighVoice')
        data.trialinfo(i,1) = 5;
    elseif strcmp(condition,'LControl')
        data.trialinfo(i,1) = 6;
    elseif strcmp(condition,'LBaseline')
        data.trialinfo(i,1) = 7;
    elseif strcmp(condition,'LLowBass')
        data.trialinfo(i,1) = 8;
    elseif strcmp(condition,'LHighVoice')
        data.trialinfo(i,1) = 9;
    end


    % Find the trials with an inconsistent number of samples
    bad_trials = find(diff(cellfun(@length, data.trial)))+1;

    % Remove the problematic trials
    % check that trialidx is valid
    if ~isempty(bad_trials)

        % create a new data structure with the specified trial removed
        % create a list of all possible indices
        all_indices = 1:length(data.trial);

        % remove the desired index from the list
        remaining_indices = setdiff(all_indices, bad_trials);

        % use the remaining indices to select the desired cells
        
        newdata = data;
        newdata.trial = data.trial(remaining_indices);
        newdata.time = data.time(remaining_indices);
        newdata.trialinfo = data.trialinfo(remaining_indices);
    else
        newdata = data;

    end

    %     % create fake sample info to help append data
    %     data.sampleinfo(i,1)=1;
    %     data.sampleinfo(i,2)=20*data.fsample;

end

