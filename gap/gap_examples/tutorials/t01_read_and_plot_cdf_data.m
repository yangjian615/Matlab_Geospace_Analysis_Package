%% Info
%   In this example we read magnetic field data from a CDF file, which was
%   provided by the Cluster Active Archive. We consider here an event
%   from 06 June 2001 at 20:00 to 21:00.
%
% Paper: 
%   M. W. Dunlop and A. Balogh. Magnetopause current as seen by Cluster. 
%   Annales Geophysicae, 23:901?907, March 2005.
%
%   The following steps are done:
%       1) Download CDF file from Cluster Science Archive
%       2) Unarchive it.
%       3) Read the file
%       4) Plot the content.
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

%% Plot data
subplot(2,1,1);             % The upper plot
plot(times, B_gse);
legend('x', 'y', 'z');      % Describes the lines on the plot
ylabel('B_{C1,GSE}');       % Describes the vertical axis
datetick('x');              % Auto-formats the x axis to display hours 
                            % and minutes

subplot(2,1,2);             % The lower plot
plot(times, pos_gse);
legend('x', 'y', 'z');
ylabel('Pos_{C1,GSE}');
xlabel('time');             % Describes the horizontal axis.
datetick('x');

