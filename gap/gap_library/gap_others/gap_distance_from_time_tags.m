function [dist_vec] = gap_distance_from_time_tags(times, pos_vec, t12)
%Computes the distance between satellite's location in two moments of time.
%USAGE:
%   [dist_vec] = gap_distance_from_time_tags(times, pos_vec, t12)
%WHERE:
%   times - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   trange - 2x1 or 1x2 vector with two time tags. Both should be given 
%           in the datenum format.
%   pos_vec - array Nx3 which has time series of position vectors.
%OUTPUT:
%   dist_vec - an array 1x3 indicating the difference vector between the
%           two positions
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Example 1:
%     t1 = datenum([2001 08 13 02 50 30]);
%     t2 = datenum([2001 08 13 02 51 00]);
%     [dist_vec] = gap_distance_from_time_tags(sel_times, c1234_pos_gse{1}, [t1, t2]);
%     vnorm(dist_vec,2)
%

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

ind = nearestpoint(t12, times);  %find array indexes corresponding to the time range
dist_vec = pos_vec(ind(2),:) - pos_vec(ind(1),:);

end
