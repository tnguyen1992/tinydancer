% -------------------------------------------------------------------------
% SUBFUNCTION which selects the appropriate artifact detection method based
% on the selected config options
% -------------------------------------------------------------------------
function [ autoart ] = artifact_detect(cfgT, data_in)

method  = cfgT.method;
sliding = cfgT.sliding;
cfgT    = removefields(cfgT, {'method', 'sliding'});

if strcmp(sliding, 'yes')                                                   % sliding window --> use own artifacts_threshold function
  autoart = artifact_sliding_threshold(cfgT, data_in);
elseif strcmp(method, 'minmax')                                             % method minmax --> use own special_minmax_threshold function
  autoart = special_minmax_threshold(cfgT, data_in);
else                                                                        % no sliding window, no minmax method --> use ft_artifacts_threshold function
  autoart = ft_artifact_threshold(cfgT, data_in);
end

end
