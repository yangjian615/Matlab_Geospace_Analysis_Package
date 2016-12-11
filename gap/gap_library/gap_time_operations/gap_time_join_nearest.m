function [time_deviation, varargout] = gap_time_join_nearest(j_times, multi_times, varargin)
%Joins time data series from various satellites or other sources. Data are 
%   "joined" only by looking for the nearest time tag.
%USAGE:
%   [time_deviation, varargout] = gap_time_join_nearest(j_times, multi_times, varargin)
%WHERE
%   j_times - is a Nx1 array with the final time tags in the datenum format (days passed 
%               from 1st January 0000). After this function, all time series will have 
%               these time tags.
%   multi_times - is a cell-array with Nsat cells, where Nsat is the number of satellites  
%               (sources), containing time tags of initial data.
%   varargin - any number of cell-arrays of data. Each cell-array must have Nsat cells
%               continaing that specific data for all satellites. Each cell must be an 
%               array Nx1 or Nx3 for scalar or vectorial data respectively. 
%OUTPUT
%   time_deviation - a cell-array with Nsat cells. Each cell is an array of the same size 
%               as j_times, indicating number of days between joined time tags and tags
%               of original sources. This is a kind of "quality check" for the time join. 
%   varargout - joined data in the same format as varargin
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%REFERENCES
%   Goetz Paschmann and Patrick W. Daly, editors. Analysis Methods for Multi-Spacecraft Data, 
%   volume SR-001 of ISSI Scientific Reports Series. ESA Publications Division, Noordwijk, 
%   The Netherlands, 1998.
%
%   Goetz Paschmann and Patrick W. Daly. Multi-Spacecraft Analysis Methods Revisited, 
%   volume SR-008 of ISSI Scientific Reports Series. ESA Publications Division, Noordwijk, 
%   The Netherlands, February 2008.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
Nsat = length(multi_times);

% Join the time axes
corr_index = cell(Nsat,1);
time_deviation = cell(Nsat,1);
for sat = 1:Nsat
    corr_index{sat} = nearestpoint(j_times, multi_times{sat});
    time_deviation{sat} = abs(multi_times{sat}(corr_index{sat}) - j_times);
end

% Join the fields according to joined time axes
nData = min(length(varargin), max(nargout,1)-1);  %-1 taken by "time_deviation"
for k = 1:nData
    multi_data = cell(size(varargin{k},1),size(varargin{k},2));
    for sat = 1:Nsat
        multi_data{sat} = varargin{k}{sat}(corr_index{sat},:);
%         multi_data{sat}(j_times < multi_times{sat}(1),:) = NaN;
%         multi_data{sat}(j_times > multi_times{sat}(end),:) = NaN;
    end
    varargout{k} = multi_data;
end

end
