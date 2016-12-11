function [sm_long, sm_lat] = gap_geo2sm_longlat(t_datenum, geo_long, geo_lat)
%Converts position from geographical coordinates (longitude, latitude)
%   to Solar Magnetic coordinates (longitude, latitude).
%USAGE:
%   [sm_long, sm_lat] = gap_geo2sm_longlat(t_datenum, geo_long, geo_lat)
%WHERE:
%   t_datenum - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   geo_long, geo_lat - geographic coordinates of the point. 
%OUTPUT:
%   sm_long, sm_lat - coordinates of the point in the Solar Magnetic coordinates. See 
%                       the gap_convert_position() function for more information about the
%                       SM coordinate system.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

xyz_geo = gap_long_lat_r2cart(geo_long, geo_lat, 1.0);  %the last one does not matter
xyz_sm = gap_convert_position(xyz_geo, t_datenum, 'GEO', 'SM');
[sm_long, sm_lat] = gap_cart2long_lat_r(xyz_sm);

end
