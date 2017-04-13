function rootPath=BOLD_LFPRootPath()
% Return the path to the root BOLD_LFP  directory
%
% This function must reside in the directory at the base of the BOLD_LFP
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(BOLD_LFP,'data')

rootPath=which('BOLD_LFPRootPath');

rootPath=fileparts(rootPath);

return
