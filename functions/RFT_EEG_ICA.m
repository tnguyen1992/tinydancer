function [ data_ica, Cmp_to_remove_all, comp] = RFT_EEG_ICA( data, cmp_nr, range, layout, data_aux, ICs_to_remove_auto )
%%
% This script runs ICA on EEG 'data' using the FT's functions.
% A number of ICs (='cmp_nr') is computed only for EEG channels (other
% channels in the dataset are ignored - although plotted). 
% 'range' and 'layout' are input for the plotting functions.
% The function works recursively - it allows you to gradually include more
% components (looking at how the data changes after they are removed). 
% Note that the components eliminated on a previous step will be
% excluded by default in all following steps (although the function performs only 1 ICA and 1 elimination of ICs). 
% Also, note that the output file includes only the EEG channels.
%
%
% FELIX: 
% 1. 'data_aux' is an ft_struct containing external electrodes in case
% you wanna plot them side by side with ICs to better interpret ICs 
% 2. 'ICs_to_remove_auto' is an array with the ICs to remove automatically
% in case you already selected them once before, and just wanna re-run the
% same pipeline without selecting them manually again
%
% You can call the functions with 'data_aux' as [] and set some
% 'ICs_to_remove_auto' if you want to removes ICs auto but not show any 
% external electrodes
% 
%% Giacomo Novembre (slightly adapted Félix)

close all

Cmp_to_remove_all   = [];
loop                = 0;

% remove non-cortical data (as indixed by 'channels' above) in order to let 'ft_rejectcomponent' work properly (otherwise you get error messages) 
cfg                 = [];
cfg.channel         = 'EEG'; % only EEG channels are considered (non EEG should not enter into ICA)
data_onlyEEG        = ft_selectdata(cfg, data);
tmp_data            = data_onlyEEG; % this is the data you can change to evaluate the effect of the ICA

% isolate non-cortical data (useful for display)
non_EEG_ch  = find(ismember(data.label,data_onlyEEG.label)==0);
cfg         = [];
cfg.channel = data.label(non_EEG_ch);
data_noEEG  = ft_selectdata(cfg, data); 

%% ICA - Run ICA, 
cfg                 = []; 
cfg.method          = 'runica';
cfg.numcomponent    = cmp_nr;        
cfg.channel         = 'EEG';
comp                = ft_componentanalysis(cfg, data_onlyEEG);

% Optimize labeling of ICs (Félix)
label = [];
for c = 1:cmp_nr    
    label{c}=['IC', num2str(c)];
end
comp.label          = label;

if exist('ICs_to_remove_auto','var')
    % Don't plot anything, only remove the ICs pre-detected
    Cmp_to_remove_all   = sort(unique(ICs_to_remove_auto));

    % remove selected components
    cfg                 = [];
    cfg.component       = Cmp_to_remove_all;
    [data_ica]          = ft_rejectcomponent(cfg, comp, data_onlyEEG);   

    display(['GIAC (felix):you have AUTOMATICALLY removed components: ' num2str(Cmp_to_remove_all) ]);
end

if ~exist('ICs_to_remove_auto','var') || isempty(ICs_to_remove_auto)

    %% 1 Plot Components
    
    % Add external electrodes if present
    if ~exist('data_aux','var') || isempty(data_aux)
        comp_PLOT = comp;
    else
        trials_aux_concat = cat(2,data_aux.trial{:});       % concatenate all trials to z-score external channels
        trials_aux_concat = trials_aux_concat - mean(trials_aux_concat,2);      % demean    
        trials_aux_concat = trials_aux_concat ./ std(trials_aux_concat,0,2);    % normalise
        trials_aux_concat = trials_aux_concat .* max(std(comp.trial{1},0,2));   % multiply by standard deviation of biggest IC (in first trial)
        [nChan , ~] = size(data_aux.trial{1,1});
        data_aux.trial = mat2cell(trials_aux_concat,nChan,cell2mat(cellfun(@(x) size(x,2),data_aux.trial,'UniformOutput',0)));   % de-concatenate back to ft_struct
        
        comp_PLOT = ft_appenddata(cfg, comp, data_aux);     % append the external channels
    end
    [ ~, ~ ] = Giac_EEGplotNanTrials( comp_PLOT, {'all'}, [], layout,  'NoReject' );
    
    %% 2 Plot Components topographies   
    figure;
    cfg                 = [];
    cfg.component       = [1:cmp_nr];       % specify the component(s) that should be plotted
    cfg.layout          = layout;
    cfg.comment         = 'no';
    ft_topoplotIC(cfg, comp);
    
    % %% 3 Plot time-locked components (to make sure you identify VW-related components)
    % cfg                 = [];
    % cfg.latency         = [-.2 +.4];
    % comp_tl             = ft_timelockanalysis(cfg, comp);
    % [ ~, ~ ] = Giac_EEGplotNanTrials( comp_tl, {'all'}, [], layout,  'NoReject' );
    
    while loop < 1
        
        %% 4 Plot Data 
        feedback_data = ft_appenddata([],tmp_data,data_noEEG);
        [ ~, ~ ]      = Giac_EEGplotNanTrials( feedback_data, {'all'}, range, layout, 'NoReject' );
        
        %% choose the components to remove
        display(['GIAC: currently you have removed components: ' num2str(Cmp_to_remove_all) ]);
        uiwait(msgbox('Press OK when you decided that some (new) components have to be removed ...', 'Rejection?'));
        Question = listdlg('PromptString','Components to remove?','SelectionMode','multiple','ListString',{'yes','no'});
    
        %% Remove components if any has been selected
        if Question == 1 % if you want to remove components ... 
    
            clear tmp_data
            
            Cmp_to_remove       = listdlg('PromptString','Choose bad components:','SelectionMode','multiple','ListString',comp.label);
            Cmp_to_remove_all   = [Cmp_to_remove_all Cmp_to_remove];
            Cmp_to_remove_all   = sort(unique(Cmp_to_remove_all));
            
            % remove selected components
            cfg                 = [];
            cfg.component       = Cmp_to_remove_all;
            [tmp_data]          = ft_rejectcomponent(cfg, comp, data_onlyEEG);   
    
            % Plot Cleaned Data 
            feedback_data = ft_appenddata([],tmp_data,data_noEEG);
            [ ~, ~ ]      = Giac_EEGplotNanTrials( feedback_data, {'all'}, range, layout, 'NoReject' );
    
            display(['GIAC: currently you have removed components: ' num2str(Cmp_to_remove_all) ]);
            uiwait(msgbox('Check your data ...', 'Rejection?'));
            Question2 = listdlg('PromptString','Select additional ICs?','SelectionMode','multiple','ListString',{'yes','no'});
            
                if Question2 == 1 % yes, other components to remove..  
                    close(4);
                    close(5);
                elseif Question2 == 2 % no, enough components removed ... 
                    loop = 1;
                end
    
        else % if no componenets have been selected, then store an empty array or overwrite an existing one
            loop = 1;
        end
    
    end % while
    
    data_ica = tmp_data;
    close all
end

end % function

