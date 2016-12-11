function [GST, slong, srasn, sdec] = gap_russell_sun_datenum(input_datenum)
%Same as russell_sun, but having inputs in datenum format (i.e. days from 1-Jan-0000)
%USAGE
%   [GST, slong, srasn, sdec] = gap_russell_sun_datenum(input_datenum)
%WHERE
%   input_datenum - an array Nx1 of datenum format (days from '1-Jan-0000') indicating
%                   the time of interest
%OUTPUT:
%   GST - Greenwich Mean Sideral Time in degrees
%   slong - ecliptic longitude
%   srasn - apparent right ascension
%   sdec - declination of the Sun in degrees
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
[year, day_of_year, sec_of_day] = gap_datenum2YandDOY(input_datenum);
[GST, slong, srasn, sdec] = gap_russell_sun(year, day_of_year, sec_of_day);

end
