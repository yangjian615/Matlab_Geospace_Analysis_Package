function [rms_diff] = gap_rms_of_difference(x, y)
%Computes Root-Mean-Square of the difference between two three-component-functions.
%USAGE:
%   [rms_diff] = gap_rms_of_difference(x, y)
%WHERE
%   x and y are arrays of row-vectors (time series of row-vectors).
%   rms_diff - if x_it is the i-th component of the vector x at time t,
%           then rms_diff = \sqrt( (\sum_t \sum_i (x-y)) / (nr time tags) )
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

rms_diff = sqrt(sum(sum((x-y).^2)) / (size(x,1)));

end

