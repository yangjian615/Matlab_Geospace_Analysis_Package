function [multi_data] = gap_time_join(j_times, time_gap, multi_times, method, input_data)
%Joins time data series from various satellites or other sources. 
%USAGE:
%   [multi_data] = gap_time_join(j_times, time_gap, multi_times, method, input_data)
%WHERE
%   j_times     - an Nx1 array with time tags in the datenum format (days passed 
%                 from 1st January 0000). 
%   time_gap    - A double specifying the maximum tolerable interval, such
%                 that T[i+1] > T[i] + time_gap is identified, and filled
%                 by linear interpolation. 
%   method      - if put to 0, a nearest-point join on gap_time_join_nearest()will be applied
%                 'boxcar': performs boxcar average on N nearest points
%                 'linear': linear interpolation between two points
%   input_data  - any number of cell-arrays of data. Each cell-array must have Nsat cells
%                 continaing that specific data for all satellites. Each cell must be an 
%                 array Nx1 or Nx3 for scalar or vectorial data respectively. 
%OUTPUT
%   multi_times - A cell-array with Nsat cells, where Nsat is the number of satellites  
%                 (sources), containing time tags of initial data.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
Nsat = length(multi_times);
multi_data = cell(size(input_data));


 % When no method is specified, apply nearest point method based on Eugen's algorithm
 if method == 0
    [time_dev, multi_data] = gap_time_join_nearest(j_times, multi_times, input_data);

 % Linear interpolation on two nearest time tags    
 elseif strcmp(method, 'linear')
     
        linear = 1;
        
        % Find the 2 nearest points to j_times[i]
        for sat = 1:Nsat
            nCol = size(input_data{sat},2); % find number of columns in each data set
            for j = 1:nCol
                for i = 1:length(j_times)
                    ind = N_nearest(j_times(i),multi_times{sat},2,linear);
                    
                    % If a time-tag in muti_times and coincide exactly with j_times
                    if length(ind)==1 
                        multi_data{sat}(i,j) = input_data{sat}(ind,j);
                        
                    else
                        dif_1 = j_times(i) - multi_times{sat}(ind(1));
                        dif_2 = j_times(i) - multi_times{sat}(ind(2));
                        
                        % If distance between one of the two nearest points and the j_times(i) exceeds time_gap threshold, 
                        % Set varargout data point = NaN
                        if (abs(dif_1) > time_gap) || (abs(dif_2) > time_gap)
                            multi_data{sat}(i,j) = NaN;
                            
                        % Else linearly interpolate between these two points
                        else
                            multi_data{sat}(i,j) = input_data{sat}(ind(1),j) + dif_1/(dif_1 - dif_2)*(input_data{sat}(ind(2),j)-input_data{sat}(ind(1),j)); 
                        end
                        
                    end
                    
                end
               
            end
        end

    
% Boxcar average of N nearest points
elseif strcmp(method, 'boxcar')
        % Determine N for Boxcar window
        %(to be worked on: find filter factor q of sample)
        % knowing that q := Nyquist freq/cutoff freq
        q = 1.35; % assumed so, for the time being;
        N = ceil(2*q);
        
        % Find the 2 nearest points to j_times[i]
        for sat = 1:Nsat
            for i = 1:length(j_times)
                
                    ind_gap = N_nearest(j_times(i),multi_times{sat}, 2, 1);
    
                    dif_1 = j_times(i) - multi_times{sat}(ind_gap(1));
                    dif_2 = j_times(i) - multi_times{sat}(ind_gap(2));
                        
                    % If distance between one of the two nearest points and the j_times(i) exceeds time_gap threshold, 
                    % Set varargout data point = NaN
                    if (abs(dif_1) > time_gap) || (abs(dif_2) > time_gap)
                        multi_data{sat}(i,:) = NaN;
                            
                    % Else linearly interpolate between these two points
                    else
                        ind = N_nearest(j_times(i),multi_times{sat}, N, 0);
                        N_nearest_data = input_data{sat}(ind,:);
                        multi_data{sat}(i,:) = mean(N_nearest_data);
                    end
                    
             end
               
            
        end
    
        
else
    % Other time re-sampling methods
end


end
