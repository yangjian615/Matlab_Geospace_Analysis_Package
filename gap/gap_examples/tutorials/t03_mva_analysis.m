%% Info
%   In this example we apply the Minimum Variance Analysis (MVA) to the
%   data from the first exercise. We know from the paper that the
%   discontinuity was observed on 06 June 2001 at 20:14
%
% Paper: 
%   M. W. Dunlop and A. Balogh. Magnetopause current as seen by Cluster. 
%   Annales Geophysicae, 23:901?907, March 2005.
%
%This code comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% Read data from CDF file:
filename = 'C1_CP_FGM_SPIN__20010611_200000_20010611_210000_V070906';
[times, B_gse, pos_gse] = gap_convert_Cluster_fgm2mat(filename);
    % times   - array Nx1 with time tags in the MatLab's datenum format
    %               (number of days since 01-January-2000). 
    % B_gse   - array Nx3 with each line having a magnetic field vector.
    % pos_gse - array Nx3 with each line having a position vector.

% Define the time interval:
tstart = datenum([2001 06 11 20 13 20]);
tend   = datenum([2001 06 11 20 15 00]);
% Select the time interval:
[sel_times, B_selected, pos_selected] ...
    = gap_select_time_range(times, [tstart tend], B_gse, pos_gse);
    % sel_times    - array Nx1 is a subset of 'times' in the specified 
    %                interval.
    % B_selected   - a subset of B_gse corresponding to time tags in
    %                sel_times.
    % pos_selected - a subset of B_gse corresponding to time tags in
    %                sel_times.

% Perform the Minimum Variance Analysis (MVA):
[n, RotMatrix, variance, B_mva] = gap_mva(B_selected);
    % n         - array 1x3 indicating the normal vector.
    % RotMatrix - rotation matrix: B_mva = B_selected * RotMatrix
    % variance  - array 3x1 with eigenvalues 
    % B_mva     - is the input field rotated into the MVA frame.

disp(n);    % Display the normal vector

%% Plot data
plot(sel_times, B_mva);
legend('min var', 'intermediate', 'max var'); 
                            % Describes the lines on the plot
ylabel('B_{C1,MVA}');       % Describes the vertical axis
xlabel('time');
datetick('x');              % Auto-formats the x axis to display hours 
                            % and minutes

                            
gap_plot_hodogram(sel_times, B_mva, [tstart tend], 'window title', variance);
% % sel_times      - time tags.
% % B_mva          - vector to be plotted.
% % [tstart tend]  - indicate the time range to be plotted.
% % 'window title' - is any information, which will appear on the plot in the
% %                  lower right corner.
% % variance       - array with eigenvalues, so that the plot can show the
% %                  ratios between variances.



