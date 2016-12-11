%% Info
%   In this example we show how to convert position vectors and magnetic 
%   field vectors from the GSE reference system to GSM.
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
    
%% Convert the reference frame
[B_gsm] = gap_convert_vector(B_gse, 0, times, 'GSE', 'GSM');
%B_gsm is a Nx3 array 
[pos_gsm] = gap_convert_position(pos_gse, times, 'GSE', 'GSM');
%pos_gsm is a Nx3 array 

%% Plot data
subplot(2,1,1);             % The upper plot
plot(times, B_gsm);
legend('x', 'y', 'z');      % Describes the lines on the plot
ylabel('B_{C1,GSM}');       % Describes the vertical axis
datetick('x');              % Auto-formats the x axis to display hours 
                            % and minutes

subplot(2,1,2);             % The lower plot
plot(times, pos_gsm);
legend('x', 'y', 'z');
ylabel('Pos_{C1,GSM}');
xlabel('time');             % Describes the horizontal axis.
datetick('x');

