function [xyz_end, rotM] = gap_convert_position(xyz_orig, times, orig_system, end_system,...
    GST, slong, srasn, sdec)
%Converts a position vector from one system to another
%USAGE:
%   [xyz_end, rotM] = gap_convert_position(xyz_orig, times, orig_system, end_system, GST, slong, srasn, sdec)  
%INPUT:   
%   xyz_orig = Nx3 position vector in the respective original system   (km)
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
%   GST, slong, srasn, sdec - optional parameters, given by gap_russell_sun_datenum(). If these
%               are not given, then they are automatically computed here.
%
% OUTPUT:
%   xyz_end = Nx3 position vector in the respective final system   (km)
%   rotM - 3x3xN contains rotation matrices for each time tag, such that
%               xyz_end(t,:) = xyz_orig(t,:) * rotM(:,:,t)
%
% See also http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
if strcmp(orig_system, end_system) == 1
    xyz_end = xyz_orig;
    rotM = repmat(eye(3),[1 1 length(times)]);
else
    MagnCollection={'GEI';'GSE';'SM';'GEO';'GSM'};
    %obtain the Sun's position
    if nargin <=4
        [GST, slong, srasn, sdec] = gap_russell_sun_datenum(times);
    end
    Nmeas = size(times,1);
    rotM = zeros(3,3,Nmeas);
    %Compare pairs of original and final systems.
    if (strcmp(orig_system, 'GEI') == 1) || (strcmp(end_system, 'GEI') == 1)
        if (strcmp(orig_system, 'GEO') == 1) || (strcmp(end_system, 'GEO') == 1)
%             A_geo2gei = [  cosd(GST), sind(GST), 0; ...
%                           -sind(GST), cosd(GST), 0; ...
%                                    0,          0, 1];
            % v_new = v_old * A_geo2gei
            %Note that the multiplication below does a transpose.
            cosd_gst = cosd(GST);
            sind_gst = sind(GST);
            rotM(3,3,:) = 1;
            rotM(1,1,:) = cosd_gst;
            rotM(2,2,:) = cosd_gst;
            if (strcmp(orig_system, 'GEO') == 1)
                rotM(1,2,:) = sind_gst;
                rotM(2,1,:) = -sind_gst;
            else
                rotM(1,2,:) = -sind_gst;
                rotM(2,1,:) = sind_gst;
            end
            xyz_end = zeros(Nmeas,3);
            for i=1:Nmeas
                xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
            end
        elseif (strcmp(orig_system, 'SM') == 1) || (strcmp(end_system, 'SM') == 1)
            %Conversion GEI->SM coordinate system
            
            %a) Find S
            s_point_gei = gap_long_lat_r2cart(srasn, sdec, 1.0);

            %b) Find D_geo (dipole direction) and convert it to D'_gei
            [Dlong_geo, Dlat_geo] = gap_geomagnetic_pole_coordinates(times);
            Dipole_geo = gap_long_lat_r2cart(Dlong_geo, Dlat_geo, 1.0);
            Dipole_gei = gap_convert_position(Dipole_geo, times, 'GEO', 'GEI', GST, slong, srasn, sdec);

            %c) Find Y and X: Y = (DxS)/|DxS|    X = YxD
            Dipole_cross_S = gap_cross(Dipole_gei, s_point_gei);
            Y = Dipole_cross_S ./ repmat(vnorm(Dipole_cross_S, 2),1,3);  %normalization
            X = gap_cross(Y, Dipole_gei);
            X = X ./ repmat(vnorm(X,2), 1, 3);  %normalization

            %d) Compute SM
            xyz_end = zeros(Nmeas,3);
            if (strcmp(orig_system, 'SM') == 1)
                for i=1:Nmeas
                    rotM(1,:,i) = X(i,:);
                    rotM(2,:,i) = Y(i,:);
                    rotM(3,:,i) = Dipole_gei(i,:);
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            else
                for i=1:Nmeas
                    rotM(:,1,i) = X(i,:)';
                    rotM(:,2,i) = Y(i,:)';
                    rotM(:,3,i) = (Dipole_gei(i,:))';
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            end
        elseif (strcmp(orig_system, 'GSM') == 1) || (strcmp(end_system, 'GSM') == 1)
            %Conversion GEI->GSM coordinate system
            
            %a) Find S
            s_point_gei = gap_long_lat_r2cart(srasn, sdec, 1.0);

            %b) Find D_geo (dipole direction) and convert it to D'_gei
            [Dlong_geo, Dlat_geo] = gap_geomagnetic_pole_coordinates(times);
            Dipole_geo = gap_long_lat_r2cart(Dlong_geo, Dlat_geo, 1.0);
            Dipole_gei = gap_convert_position(Dipole_geo, times, 'GEO', 'GEI', GST, slong, srasn, sdec);

            %c) Find Y and X: Y = (DxS)/|DxS|    X = YxD
            Dipole_cross_S = gap_cross(Dipole_gei, s_point_gei);
            Y = Dipole_cross_S ./ repmat(vnorm(Dipole_cross_S, 2),1,3);  %normalization
            Z = gap_cross(s_point_gei, Y);
            Z = Z ./ repmat(vnorm(Z,2), 1, 3);  %normalization

            %d) Compute SM
            xyz_end = zeros(Nmeas,3);
            if (strcmp(orig_system, 'GSM') == 1)
                for i=1:Nmeas
                    rotM(1,:,i) = s_point_gei(i,:);
                    rotM(2,:,i) = Y(i,:);
                    rotM(3,:,i) = Z(i,:);
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            else
                for i=1:Nmeas
                    rotM(:,1,i) = (s_point_gei(i,:))';
                    rotM(:,2,i) = (Y(i,:))';
                    rotM(:,3,i) = (Z(i,:))';
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            end
        elseif (strcmp(orig_system, 'GSE') == 1) || (strcmp(end_system, 'GSE') == 1)
            %Conversion GEI->GSE coordinate system
            %a) Find S
            s_point_gei = gap_long_lat_r2cart(srasn, sdec, 1.0);

            %b) Find the ecliptic pole
            ecliptic_pole_gei = [0, -0.3981, 0.9173];
            ecliptic_pole_gei = ecliptic_pole_gei / norm(ecliptic_pole_gei);

            %c) Find Y = ecliptic_pole x S
            Y = gap_cross(repmat(ecliptic_pole_gei, Nmeas, 1), s_point_gei);
            Y = Y ./ repmat(vnorm(Y,2), 1, 3);  %normalization

            %  correct the ecliptic pole, because it is not perfect perpendicular to s_point
            ecliptic_pole_gei = gap_cross(s_point_gei, Y);
            ecliptic_pole_gei = ecliptic_pole_gei ./ repmat(vnorm(ecliptic_pole_gei,2), 1, 3);  %normalization
            
            %d) Compute SM
            xyz_end = zeros(Nmeas,3);
            if (strcmp(orig_system, 'GSE') == 1)
                for i=1:Nmeas
                    rotM(1,:,i) = s_point_gei(i,:);
                    rotM(2,:,i) = Y(i,:);
                    rotM(3,:,i) = ecliptic_pole_gei(i,:);
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            else
                %ecliptic_pole_gei = ecliptic_pole_gei';
                for i=1:Nmeas
                    rotM(:,1,i) = (s_point_gei(i,:))';
                    rotM(:,2,i) = (Y(i,:))';
                    rotM(:,3,i) = (ecliptic_pole_gei(i,:))';
                    xyz_end(i,:) = xyz_orig(i,:) * rotM(:,:,i);
                end
            end
        end
    elseif (IsOneFromCollection(orig_system, MagnCollection)) && (IsOneFromCollection(end_system, MagnCollection))
        [xyz_gei, rotM1] = gap_convert_position(xyz_orig, times, ...
            orig_system, 'GEI', GST, slong, srasn, sdec);
        [xyz_end, rotM2] = gap_convert_position(xyz_gei, times, ...
            'GEI', end_system, GST, slong, srasn, sdec);
        for i=1:Nmeas
            rotM(:,:,i) = rotM1(:,:,i) * rotM2(:,:,i);
        end
    else
        error('This pair of reference systems are not yet supported.');
    end
end

end  %main function end

%=========================================================================================
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

