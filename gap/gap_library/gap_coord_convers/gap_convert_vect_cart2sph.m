function [vec_end] = gap_convert_vect_cart2sph(vec_orig, longlat)
%Converts a vector from the cartesian reference 
%   system to a spherical local reference system (North, East, Center). NEC
%   frame is a frame whose orientation depends on the location of the point
%   of interest. NEC means positive to North, positive to East, positive
%   toward center. The final cartesian frame is universal for all points on
%   the sphere, has the origin at the center of the sphere and axes as
%   'z' axis is directed directly to North (latitude 90 degrees)
%   'x' axis is in direction latitude and longitude=0 degrees
%   'y' completes the triad.
%USAGE:
%   [vec_end] = gap_convert_vect_cart2sph(vec_orig, longlat)
% 
%INPUT:   
%   vec_orig = Nx3 an array of N vector in the respective original cartesian system.
%   longlat = Nx2 array of N pair of longitude and latitudes IN DEGREES.
% OUTPUT:
%   vec_end = Nx3 vector in the NEC frame.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
Nmeas = size(vec_orig,1);
col = 90.0 - longlat(:,2);  %colatitude
T = zeros(3,3,Nmeas);
%T = [sind(col)*cosd(long), sind(col)*sind(long),  cosd(col); ...
%    cosd(col)*cosd(long), cosd(col)*sind(long), -sind(col); ...
%    -sind(long),           cosd(long),         0];
sind_col = sind(col);
cosd_col = cosd(col);
sind_long = sind(longlat(:,1));
cosd_long = cosd(longlat(:,1));
T(1,1,:) = sind_col .* cosd_long;
T(2,1,:) = sind_col .* sind_long;
T(3,1,:) = cosd_col;
T(1,2,:) = cosd_col .* cosd_long;
T(2,2,:) = cosd_col .* sind_long;
T(3,2,:) = -sind_col;
T(1,3,:) = -sind_long;
T(2,3,:) = cosd_long;
T(3,3,:) = 0;
vec_end = zeros(Nmeas,3);
vec_spheric = zeros(Nmeas,3);
for i=1:Nmeas
    vec_spheric(i,:) = vec_orig(i,:) * T(:,:,i);
end
vec_end(:,1) = -vec_spheric(:,2);  % colatitude is directed to South. nec(1) is to North
vec_end(:,2) =  vec_spheric(:,3);  % East.
vec_end(:,3) = -vec_spheric(:,1); % toward center; R is away from center.

end
