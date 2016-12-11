function [times, B_gse_vec, pos_gse_vec] = gap_convert_Cluster_fgm2mat(cdf_filename, mat_filename)
%Reads a CDF file with FGM data from Cluster Active Archive.
%Usage:
%   [times, B_gse_vec, pos_gse_vec] = gap_convert_Cluster_fgm2mat(cdf_filename, mat_filename)
%WHERE
%   cdf_filename - is a string with the path of the FGM data in a CDF file.
%   mat_filename - (optional) the name of the .mat file, where the content
%               will be saved.
%OUTPUT
%   times     - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   B_gse_vec - is a Nx3 array with rows-vectors representing magnetic
%               field in the coordinates provided by CDF file.
%   pos_gse_vec - is a Nx3 array with rows-vectors representing position
%               vectors of the satellite, in the coordinates provided by CDF file.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
data = cdfread(cdf_filename, 'ConvertEpoch', true);
times = cell2mat(data(:,1));
B_gse_vec = cell2mat(data(:,3)')';
pos_gse_vec = cell2mat(data(:,5)')';
if (nargin > 1) && ischar(mat_filename) && (~isempty(mat_filename))
    save(mat_filename, 'times', 'B_gse_vec', 'pos_gse_vec');
end

end
