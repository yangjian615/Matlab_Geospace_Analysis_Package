function [GST, slong, srasn, sdec] = gap_russell_sun(year, day_year, sec_day)
%The calculation of the position of the sun
%USAGE:
%   [GST, slong, srasn, sdec] = gap_russell_sun(year, day_year, sec_day)
%WHERE
%   year - an array Nx1 of integers representing the actual year in the gregorian calendar (e.g. 2010).
%   day_year - an array Nx1 of integers 1..356
%   sec_day - an array Nx1 of integers 0..86400
%OUTPUT:
%   GST - Greenwich Mean Sideral Time in degrees
%   slong - ecliptic longitude
%   srasn - apparent right ascension
%   sdec - declination of the Sun in degrees
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% REFERENCE
% Taken from Russell. Original publication:
% Cosmic Electrodynamics, 2, 184-196, 1971. All rights reserved. 
% Copyright 1971 by D. Reidel, Publishing Company Dordrecht-Holland
% http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
rad = 57.29578;
fday = sec_day / 86400;
DJ = 365 * (year - 1900) + (year - 1901)/4 + day_year + fday - 0.5;
t = DJ / 36525;
vl =  gap_remainder(279.696678 + 0.9856473354*DJ,                        360);
GST = gap_remainder(279.690983 + 0.9856473354 * DJ + 360.0*fday + 180.0, 360);
g =   gap_remainder(358.475845 + 0.985600267*DJ,                         360) / rad;
slong = vl + (1.91946 - 0.004789 * t) .* sin(g) + 0.020094 * sin(2.0*g);
obliq = (23.45229 - 0.0130125 * t) / rad;
slp = (slong - 0.005686) / rad;
sindd = sin(obliq) .* sin(slp);
cosdd = sqrt(1.0 - sindd .^2);
sdec = rad * atan(sindd ./ cosdd);
srasn = 180.0 - rad * atan2(sindd ./ cosdd ./ tan(obliq), -cos(slp) ./ cosdd);

end
