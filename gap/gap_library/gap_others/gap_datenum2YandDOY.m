function [year, day_of_year, sec_of_day] = gap_datenum2YandDOY(serial_date_number)
%Extracts year, day of the year and seconds of the day from the matlab's 
%   native serial date number.
%USAGE: 
%   [year, day_of_year, sec_of_day] = datenum2YandDOY(serial_date_number)
%WHERE
%   serial_date_number -- a Nx1 array of datenum formats (days passed from 1st January 0000)
%OUTPUT
%   year -- a Nx1 array of integers (real Gregorian years)
%   day_of_year -- a Nx1 array of integers
%   sec_of_day -- a Nx1 array of integers
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

time_vec = datevec(serial_date_number);
Nentries = size(time_vec,1);
timevec_yymmdd_only = time_vec;
timevec_yymmdd_only(:,4:6) = zeros(Nentries,3);  %make hours, minutes, seconds = 0
timevec_year_only = timevec_yymmdd_only;
timevec_year_only(:,2) = zeros(Nentries,1);  %make month = 0
timevec_year_only(:,3) = ones(Nentries,1);  %days must be 1
year = time_vec(:,1);
datenum_yymmdd_only = datenum(timevec_yymmdd_only);
day_of_year = datenum_yymmdd_only - datenum(timevec_year_only);
sec_of_day = (serial_date_number - datenum_yymmdd_only) * 86400;  %convert to seconds

end
