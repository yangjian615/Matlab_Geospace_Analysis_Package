function [B_end_frame] = gap_igrf(gh, times, xyz, xyz_frame, end_frame, itype)
%Gives the IGRF magnetic field model in specified coordinates
%  USAGE:   [B]=igrf11syn(fyears,alt,nlat,elong) 
% 
%  INPUT:   gh - model's coefficients as given by the function gap_prepare_igrf()
%           times - Nx1 date time in datenum format (days since 01 January 0000) 
%           xyz    - Nx3 position vector in the respective system (end_frame) (km)
%           xyz_frame - a string indicate the frame where xyz is defined. Can be one of
%                   the following: 'GEI';'GSE';'SM';'GEO';'GSM';'NEC'.
%                   If this is set to NEC, then xyz is actually a vector of type: [long, lat, r].
%                   See the function gap_convert_position() for more information.
%           end_frame - a string indicate the frame to which B field must be converted.
%                   Can be one of the following: 'GEI';'GSE';'SM';'GEO';'GSM';'NEC'
%                   See the function gap_convert_vector() for more information.
%           itype  - is considered only in NEC system. Otherwise 2 by default.
%                  = 1 if geodetic (spheroid)
%                  = 2 if geocentric (sphere)
% OUTPUT:
%           B      = Nx3 vector in the respective coordinate system (nano Teslas)
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
error(nargchk(5,6,nargin)) ;  %only itype is optional
%Default values:
if nargin <6
    itype = 2;  %default is from the center of the earth.
elseif (itype == 1) && (strcmp(xyz_frame, 'NEC') ~= 1)
    error('itype=1 does not make sense for a frame different from NEC.');
end

N_meas = size(times, 1);
if N_meas ~= size(xyz,1)
    error('Time and position vectors must be of the same length');
end
if (size(times, 2)~=1) || (size(xyz,2)~=3)
    error('fyears must be a Nx1 array. xyz must be a Nx3 array.');
end

if (strcmp(xyz_frame, 'NEC') == 1)
    long = xyz(:,1);
    lat = xyz(:,2);
    radius = xyz(:,3);
else
    [xyz_geo] = gap_convert_position(xyz, times, xyz_frame, 'GEO');

    %IGRF works with geographic latitudes and longitudes. We also need this to convert to further frames.
    [long, lat, radius] = gap_cart2long_lat_r(xyz_geo);
end

%Get IGRF data in NEC frame
[year, day_of_year] = gap_datenum2YandDOY(times);
fyears = year + (day_of_year / 365.25);
B_nec = zeros(N_meas, 3);
parfor i=1:N_meas
    B_nec(i,:) = crino_igrf11syn(gh,fyears(i),itype,radius(i),lat(i),long(i));
end

%Convert the result to the desired reference system
[B_end_frame] = gap_convert_vector(B_nec, [long, lat], times, 'NEC', end_frame);

end

