function [selected_times, varargout] = gap_select_time_range(times, trange, varargin)
%Selects a subset of data based on the desired time range.
%USAGE:
%   [selected_times, varargout] = gap_select_time_range(times, trange, varargin)
%WHERE:
%   times - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   trange - 2x1 or 1x2 vector with two elements: start and end of the time span. Both
%               should be given in the datenum format.
%   variable input - any number of data variables in the format NxM, so that the rows correspond
%               to time tags given in the vector 'times'. The code does not check this.
%OUTPUT:
%   variable output - data variables which correspond to the selected time tags.
%EXAMPLE 1:
%   [stimes] = gap_select_time_range((1:10)', [4 6])
%   RESULT:
%         stimes =
% 
%              4
%              5
%              6
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%EXAMPLE 2:
%   times = (1:10)';
%   [stimes, A, B] = gap_select_time_range(times, [4;6], [times*2, sind(times)], cosd(times))
%   RESULT:
%         stimes =
% 
%              4
%              5
%              6
% 
% 
%         A =
% 
%             8.0000    0.0698
%            10.0000    0.0872
%            12.0000    0.1045
% 
% 
%         B =
% 
%             0.9976
%             0.9962
%             0.9945

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
ind = nearestpoint(trange, times);  %find array indexes corresponding to the time range
selected_times = times(ind(1):ind(2));  %copy time tags in the time range

% nVarargsIn = length(varargin);
% nVarargsOut = max(nargout,1) - 1;  %1 out argument is taken by selected_times (taken from
                %examples from Matlab's help.
  
nData = min(length(varargin), max(nargout,1) - 1);
for k = 1:nData
    varargout{k} = varargin{k}(ind(1):ind(2),:);  %all data 
end               

end
