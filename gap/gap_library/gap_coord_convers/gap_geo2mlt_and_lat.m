function [mlt, sm_lat] = gap_geo2mlt_and_lat(t_datenum, geo_long, geo_lat)
%Converts position from geographical coordinates (longitude, latitude)
%   to magnetic local time (MLT) and Solar Magnetic (SM) latitude.
%USAGE:
%   [mlt, sm_lat] = gap_geo2mlt_and_lat(t_datenum, geo_long, geo_lat)
%WHERE:
%   t_datenum - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   geo_long, geo_lat - geographic coordinates of the point. 
%OUTPUT:
%   mlt    - magnetic local time. It is equal to midnight when sm_long=180 and is noon at sm_long=0.
%   sm_lat - coordinates of the point in the Solar Magnetic coordinates. See the 
%                       gap_convert_position() function for more information about the
%                       SM coordinate system.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

%First part is duplicate of the function 
%[sm_long, sm_lat] = gap_geo2sm_longlat(t_datenum, geo_long, geo_lat)
xyz_geo = gap_long_lat_r2cart(geo_long, geo_lat, 1.0);  %the last one does not matter
xyz_sm = gap_convert_position(xyz_geo, t_datenum, 'GEO', 'SM');
[sm_long, sm_lat] = gap_cart2long_lat_r(xyz_sm);

%Do a correction of 12 hours. SM has X axix pointing to the sun.
%   We need however, 0 degrees to be exactly opposite to the sun.
mlt = gap_remainder(sm_long + 180.0, 360.0) * 24.0 / 360.0;
end
