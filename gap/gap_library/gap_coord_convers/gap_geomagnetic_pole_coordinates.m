function [geo_long, geo_lat] = gap_geomagnetic_pole_coordinates(t_datenum)
%Computes the coordinates of geomagnetic north pole from IGRF data from many years. 
%USAGE:
%   [geo_long, geo_lat] = gap_geomagnetic_pole_coordinates(t_datenum)
%WHERE:
%   t_datenum - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%OUTPUT:
%   geo_long, geo_lat - geographic coordinates of the North geomagnetic pole (different from
%                       the N geographic or the N magnetic pole). Ranges:
%                       geo_long - 0..360 from GMT eastward
%                       geo_lat - -90..90 with positive northward.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

% IGRF_N_poles = ...
%     [1970 ,  1975,  1980,  1985,  1990,  1995,  2000,  2005,  2010,  2015; ...
%      78.66, 78.76, 78.88, 79.04, 79.21, 79.39, 79.61, 79.82, 80.08, 80.36; ...
%     -70.18,-70.47,-70.76,-70.90,-71.13,-71.42,-71.57,-71.81,-72.22,-72.62];

%Next is in datenum format. (time;geo lat;geo long)
IGRF_N_poles = ...
    [719529, 721355, 723181, 725008, 726834, 728660, 730486, 732313, 734139, 735965; ...
     78.66, 78.76, 78.88, 79.04, 79.21, 79.39, 79.61, 79.82, 80.08, 80.36; ...
    -70.18,-70.47,-70.76,-70.90,-71.13,-71.42,-71.57,-71.81,-72.22,-72.62]';

%temp_years = t_datenum / 365.23;  % A kind of datevec to extract years only
if (~isempty(find((t_datenum >= max(IGRF_N_poles(:,1))) | (t_datenum <= min(IGRF_N_poles(:,1))), 1)))
    error('Only year between 1970 and 2015 are supported!');
end

%xmx1 == x minus x1 = x-x1
[ind,xmx1] = nearestpoint(t_datenum, IGRF_N_poles(:,1), 'previous');
Dlatlong_y1 = IGRF_N_poles(ind,2:3);

%x2mx == x2 minus x = x2-x
[ind,x2mx] = nearestpoint(t_datenum, IGRF_N_poles(:,1), 'next');
Dlatlong_y2 = IGRF_N_poles(ind,2:3);

% y = y1 + (y2 - y1) * [(x - x1) / (x2 - x1)] where  t = [(x - x1) / (x2 - x1)]
t = xmx1 ./ (x2mx + xmx1);
Dlatlong = Dlatlong_y1 + (Dlatlong_y2 - Dlatlong_y1) .* repmat(t,1,2);
geo_lat = Dlatlong(:,1);
geo_long = Dlatlong(:,2);
ind = find(geo_long < 0);  %make all degrees in the range 0..360
geo_long(ind) = geo_long(ind) + 360;  

end
