function [xyz] = gap_long_lat_r2cart(long, lat, radius)
%Converts longitude, latitude and radius to xyz type. Analogy to sph2cart().
%USAGE:
%   [xyz] = gap_long_lat_r2cart(long, lat, radius)
%WHERE
%   long, lat - longitude and latitude given in degrees (can be in any frame). Each must be
%               an array Nx1.
%   radius - an array Nx1 of distances from the origin of coordinates
%OUTPUT
%   xyz - is a Nx3 array with data of cartesian type
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
xyz = zeros(length(long),3);
cosd_latitude = radius .* cosd(lat);
xyz(:,1) = cosd_latitude .* cosd(long);
xyz(:,2) = cosd_latitude .* sind(long);
xyz(:,3) = radius .* sind(lat);

end
