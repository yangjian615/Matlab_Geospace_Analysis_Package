function gh_coeffs = gap_prepare_igrf()
%gap_prepare_igrf   Loads coefficients needed to compute the IGRF'11 model.
%
%Syntax
%   gh_coeffs = gap_prepare_igrf();
%
%Format:
%   An 1x3255 array where coefficients from each year (155 coeffs/year) are placed after
%   the previous year. It starts from 1900 and goes each 5 years until 2010. The last are
%   secular variations for 2010 - 2015.
%
%See also:
%   http://www.mathworks.com/matlabcentral/fileexchange/28874-igrf-magnetic-field
%   http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

% load gap_igrf_gh_coefficients.mat;
% gh_coeffs = gap_igrf_GH_coefficients;
% clear gap_igrf_GH_coefficients;

load ../other_toolboxes/IGRF11/IGRF/GHcoefficients.mat;
gh_coeffs = gh;
clear gh;

end

