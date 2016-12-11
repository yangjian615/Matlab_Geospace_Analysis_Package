function [long, lat, radius] = gap_cart2long_lat_r(xyz)
%Converts xyz type to longitude, latitude and radius. Analogy to cart2sph().
%USAGE:
%   [long, lat, radius] = gap_cart2long_lat_r(xyz)
%WHERE
%   xyz - is a Nx3 array with data of cartesian type
%OUTPUT
%   long   - array Nx1 of longitudes given in degrees (0..360).
%   lat    - array Nx1 of latitudes given in degrees (-90..90 with positive to geographic North).
%   radius - an array Nx1 of distances from the origin of coordinates
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

[long, lat, radius] = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));
lat = lat * 180.0 / pi;
long = long * 180.0 / pi;
ind = find(long < 0);  %make all degrees in the range 0..360
long(ind) = long(ind) + 360;

end
