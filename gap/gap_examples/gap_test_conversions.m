function  gap_test_conversions()
%GAP_TEST_CONVERSIONS Tests the gap_convert_position() function
%OUTPUT:
% Difference SM from GSE: 0.003960
% Difference GEI from SM: 0.004575
% Difference GSM from SM: 0.002665
% Difference GSE end - original: 0.000000
%NOTE that the last should be 0 because it is return to the initial state. If it is not 0,
%   then the determinant of the rotation matrix is not perfect 1.0
%
%This code comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

times = [7.3480e+05; 7.3480e+05 - 365.25; 7.3480e+05 - 60; 7.3480e+05 - 3652.5];
xyz_gse =  [12000.0 , 0.0       , 0.0;      ...
            0.0     , 12000.0   , 0.0;      ...
            0.0     , 0.0       , 12000.0;  ...
            7000.0  , 2000.0    , 8000.0];

[xyz_sm] = gap_convert_position(xyz_gse, times, 'GSE', 'SM');
xyz_sm_0 =  1.0e+04 * ...
  [ 1.1601    0.0000   -0.3069; ...
   -0.0128    1.1471   -0.3521; ...
   -0.1659    0.2350    1.1650; ...
    0.7677    0.5770    0.4976];
dif = evaluate_diff_between_matrices(xyz_sm, xyz_sm_0);
s = sprintf('Difference SM from GSE: %f', dif);
disp(s);
    
[xyz_gei] = gap_convert_position(xyz_sm, times, 'SM', 'GEI');
xyz_gei_0 =  1.0e+04 * ...
  [-1.0421   -0.5460   -0.2367; ...
    0.5904   -0.9584   -0.4159; ...
   -0.0002   -0.4776    1.1008; ...
   -0.5138   -0.7915    0.5287];
dif = evaluate_diff_between_matrices(xyz_gei, xyz_gei_0);
s = sprintf('Difference GEI from SM: %f', dif);
disp(s);

[xyz_gsm] = gap_convert_position(xyz_sm, times, 'SM', 'GSM');
xyz_gsm_0 =  1.0e+04 * ...
  [ 1.2000    0.0000    0.0000; ...
   -0.0000    1.1471   -0.3523; ...
    0.0000    0.2350    1.1768; ...
    0.7000    0.5770    0.5891];
dif = evaluate_diff_between_matrices(xyz_gsm, xyz_gsm_0);
s = sprintf('Difference GSM from SM: %f', dif);
disp(s);

[xyz_gse_conv] = gap_convert_position(xyz_gsm, times, 'GSM', 'GSE');
dif = evaluate_diff_between_matrices(xyz_gse_conv, xyz_gse);
s = sprintf('Difference GSE end - original: %f', dif);
disp(s);

end

%========================================================================================

function result = evaluate_diff_between_matrices(A, B)

dif = abs(A - B);
dif = sqrt(sum(sum(dif .^ 2)) / (size(A,1) * size(A,2)));
Anorm = sqrt(sum(sum(A .^ 2)) / (size(A,1) * size(A,2)));
result = (dif * 100.0) / Anorm;
end
