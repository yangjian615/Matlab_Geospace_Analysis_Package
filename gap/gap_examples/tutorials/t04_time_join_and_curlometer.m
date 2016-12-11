%% Info
%   In this example we show how to join the data and compute the current
%   from four satellites using the curlometer method.
%
% Paper of the event: 
%   M. W. Dunlop and A. Balogh. Magnetopause current as seen by Cluster. 
%   Annales Geophysicae, 23:901?907, March 2005.
%
%This code comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.


%% Read data from CDF files:
filename1 = 'C1_CP_FGM_SPIN__20010611_200000_20010611_210000_V070906.cdf';
filename2 = 'C2_CP_FGM_SPIN__20010611_200000_20010611_210000_V070906.cdf';
filename3 = 'C3_CP_FGM_SPIN__20010611_200000_20010611_210000_V070906.cdf';
filename4 = 'C4_CP_FGM_SPIN__20010611_200000_20010611_210000_V070906.cdf';
[times1, B1, pos1] = gap_convert_Cluster_fgm2mat(filename1);
[times2, B2, pos2] = gap_convert_Cluster_fgm2mat(filename2);
[times3, B3, pos3] = gap_convert_Cluster_fgm2mat(filename3);
[times4, B4, pos4] = gap_convert_Cluster_fgm2mat(filename4);

% Group variables into cells
c1234_times = {times1, times2, times3, times4};
c1234_B = {B1, B2, B3, B4};
c1234_pos = {pos1, pos2, pos3, pos4};


%% Join time series of data
joined_times = c1234_times{1};
time_gap = 10/24.0/3600.0; % tolerate maximum 10 seconds time gap
c1234joined_pos = ...
    gap_time_join(joined_times, time_gap, c1234_times, 'linear', c1234_pos);
c1234joined_B = ...
    gap_time_join(joined_times, time_gap, c1234_times, 'linear', c1234_B);


%% Currents: Curlometer
[curlB, K_recipr, divB] = gap_curlometer(c1234joined_pos, c1234joined_B);
J = curlB * 1e-12 / (1.25663706e-06);  %Convert curl into currents,
                                       % from nT/km to A/m^2.
J = J * 1e9;    % Convert from A/m^2 to nA/m^2.


%% Plot data
subplot(3,1,1);             % The upper plot
plot(joined_times, J);
legend('x', 'y', 'z');
ylim([-10 10]);             % Some limits, taken from the paper.
ylabel('J_{GSE} 10^{-9} Am^{-2}');
xlabel('time');             % Describes the horizontal axis.
datetick('x');              % Auto-formats the x axis to display hours 
                            % and minutes

subplot(3,1,2);             % The second plot
plot(joined_times, abs(divB ./ vnorm(curlB,2)) * 100.0);
                            % vnorm(B,2) computes norms of vectors B
ylim([0 200]);              % Some limits, taken from the paper.
ylabel('divB %');
datetick('x');

subplot(3,1,3);             % The third plot
plot(times1, vnorm(B1,2), times2, vnorm(B2,2) ...
    , times3, vnorm(B3,2), times4, vnorm(B4,2));
legend('C1', 'C2', 'C3', 'C4');  % label the lines
ylim([0 30]);               % Some limits, taken to look good. 
ylabel('|B| nT');       
datetick('x');     


