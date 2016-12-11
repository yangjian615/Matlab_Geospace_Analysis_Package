function [vec_end] = gap_convert_vector(vec_orig, longlat, times, orig_system, end_system,...
    GST, slong, srasn, sdec)
%Converts a vector from one reference system to another.
%USAGE:
%   [vec_end] = gap_convert_vector(vec_orig, longlat, times, orig_system, end_system, GST, slong, srasn, sdec) 
% 
%INPUT:   
%   vec_orig = Nx3 an array of N vector in the respective original system.
%   longlat = Nx2 array of N pair of longitude and latitudes. This variable makes sense only
%               only if either orig_system or end_system is in 'NEC' frame. Otherwise this 
%               can be any value.
%   times = Nx1 array with time values (as in datenum format, i.e. fractional days
%                    from 'Jan-1-0000 00:00:00'). This is needed for most conversions
%                    because most ref systems depend on the relative position of the Sun.
%                    and of the rotation axis. 
%   orig_system, end_system = strings describing a reference systems. Currently
%               the following systems are supported:
%               'GSE'  - Geocentric Solar Ecliptic: has its X-axis pointing from the Earth
%                           towards the Sun and its Y-axis is chosen to be in the ecliptic
%                           plane pointing towards dusk (thus opposing planetary motion). 
%                           Its Z-axis is parallel to the ecliptic pole. 
%               'GEO'  - Geographic (maybe same as ITRF): its X-axis is in the Earth's 
%                           equatorial plane but is fixed with the rotation of the Earth 
%                           so that it passes through the Greenwich meridian (0o longitude).
%                           Its Z-axis is parallel to the rotation axis of the Earth, and 
%                           its Y-axis completes a right-handed orthogonal set (Y= ZxX). 
%               'GEI'  - Geocentric Equatorial Inertial: has its X-axis pointing from the 
%                           Earth towards the first point of Aries (the position of the 
%                           Sun at the vernal equinox). This direction is the intersection
%                           of the Earth's equatorial plane and the ecliptic plane and 
%                           thus the X-axis lies in both planes. The Z-axis is parallel 
%                           to the rotation axis of the Earth and Y completes the 
%                           right-handed orthogonal set.
%               'SM'   - Solar Magnetic: Z-axis is chosen parallel to the north magnetic 
%                           pole and the Y-axis perpendicular to the Earth-Sun line 
%                           towards dusk. note that in this system the X-axis does not 
%                           point directly at the Sun.
%               'GSM'  - Geocentric Solar Magnetospheric: has its X-axis from the Earth 
%                           to the Sun. The Y-axis is defined to be perpendicular to the 
%                           Earth's magnetic dipole so that the X-Z plane contains the 
%                           dipole axis. The positive Z- axis is chosen to be in the same 
%                           sense as the northern magnetic pole. 
%               'NEC'  - (Northward, Eastward, Centroid) frame is a local frame with the origin in 
%                           the geometric center of the vector feedback magnetometer 
%                           (VFM). The radial component points from the centre of the VFM 
%                           towards the centre of the Earth (defined in ITRF). The 
%                           North (N) and East (E) components point from the center of the 
%                           VFM towards North and East, i.e. along the local tangent to the 
%                           meridian, respectively the parallel, of the sphere (defined in 
%                           ITRF) with radius 1from the center of the Earth to the center 
%                           of the instrument.
%   GST, slong, srasn, sdec - optional parameters, given by gap_russell_sun_datenum(). If these
%               are not given, then they are automatically computed here.
%
% OUTPUT:
%   vec_end = Nx3 position vector in the respective final system   (km)
%
% See also http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html
%   and the footnotes of the page 31 of the pdf file: esamultimedia.esa.int/docs/EEUCM/SWARM_TPA.pdf
%   for the definition of NEC frame.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
if strcmp(orig_system, end_system) == 1
    vec_end = vec_orig;
    %rotM = repmat(eye(3),[1 1 length(times)]);
else
    MagnCollection={'GEI';'GSE';'SM';'GEO';'GSM'};
    %obtain the Sun's position
    if nargin <=5
        [GST, slong, srasn, sdec] = gap_russell_sun_datenum(times);
    end
    Nmeas = size(times,1);
    %rotM = zeros(3,3,Nmeas);
    %Compare pairs of original and final systems.
    if (IsOneFromCollection(orig_system, MagnCollection)) && (IsOneFromCollection(end_system, MagnCollection))
        %xyz_orig = zeros(Nmeas, 3); -- we don't need to store this variable
        [xyz_end, rotM] = gap_convert_position(zeros(Nmeas, 3), times, orig_system, end_system,...
            GST, slong, srasn, sdec);
        vec_end = zeros(Nmeas,3);
        for i=1:Nmeas
            vec_end(i,:) = vec_orig(i,:) * rotM(:,:,i);
        end
    elseif (strcmp(orig_system, 'GEO') == 1) && (strcmp(end_system, 'NEC') == 1)
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
    elseif (IsOneFromCollection(orig_system, MagnCollection)) && (strcmp(end_system, 'NEC') == 1)
        [vec_geo] = gap_convert_vector(vec_orig, longlat, times, orig_system, 'GEO',...
            GST, slong, srasn, sdec);
        [vec_end] = gap_convert_vector(vec_geo, longlat, times, 'GEO', end_system,...
            GST, slong, srasn, sdec);
    elseif (strcmp(end_system, 'GEO') == 1) && (strcmp(orig_system, 'NEC') == 1)
        vec_spheric = zeros(Nmeas,3);
        vec_spheric(:,2) = -vec_orig(:,1);
        vec_spheric(:,3) =  vec_orig(:,2);
        vec_spheric(:,1) = -vec_orig(:,3);

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
        for i=1:Nmeas
            vec_end(i,:) = vec_spheric(i,:) * (T(:,:,i))';  %an easy fix, but computationally expensive
                                                       % to do transpose every time.
        end
    elseif (IsOneFromCollection(end_system, MagnCollection)) && (strcmp(orig_system, 'NEC') == 1)
        [vec_geo] = gap_convert_vector(vec_orig, longlat, times, orig_system, 'GEO',...
            GST, slong, srasn, sdec);  %NEC -> GEO
        [vec_end] = gap_convert_vector(vec_geo, longlat, times, 'GEO', end_system,...
            GST, slong, srasn, sdec);  %GEO -> wanted system
    else
        error('This pair of reference systems are not yet supported.');
    end
end

end  %main function end

%% =========================================================================================
function [is_in_collection] = IsOneFromCollection(str_system, collection)
%IsOneFromCollection Checks if a string is included in a set of strings
%USAGE:
%   [is_in_collection] = IsOneFromCollection(str_system, collection)
% 
%INPUT:
%   str_system -- a string which must be looked for in the collection
%   collection -- an array of cells of type strings (e.g. A = {'gse'; 'gsm'; 'sm'}
%OUTPUT:
%   is_in_collection -- bool is true if str_system is part of collection

is_in_collection = ~isempty(find(strcmp(str_system, collection)>0, 1));

end
