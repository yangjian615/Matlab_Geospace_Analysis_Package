function [normal, RotMatrix, variance, B_mva, nr_points, Delta_phi] ...
    = gap_mva_with_time_range(times, time_range, B)
%Performs the unconstrained Minimum Variance Analysis (MVA) on the subset of the
%       time series B, based on the time range
%USAGE:
%   [n, RotMatrix, variance, B_mva] = gap_mva_with_time_range(times, time_range, B)
%WHERE
%   times - Nx1 array with times in the datenum format (days since "1-Jan-0000").
%   time_range - 2x1 or 1x2 vector with two elements: start and end of the time span. Both
%           should be given in the datenum format. Only this interval is used in MVA.
%   B - is a Nx3 array with rows as cartesian vectors
%OUTPUT
%   normal   - vector 1x3 representing the normal vector in the original frame.
%   RotMatrix - 3x3 matrix representing the rotation matrix between the original and the MVA
%           frame. The conversion is done as follows:
%           B_mva = B * RotMatrix
%           provided that both B_mva and B are Nx3 time series of vectors.
%           The maximum variance direction is RotMatrix(:,3)'
%           The normal direction is RotMatrix(:,1)'
%   variance - 3x1 vector with the tree eigenvalues of the MVA problem. These eigenvalues
%           represent variance of the field in the x, y and z direction in the MVA frame.
%           variance(1) is the smallest variance. The x coordinate in the MVA has the smallest
%           change.
%   B_mva - optional array Nx3 with ALL vectors B transformed into the MVA frame. 
%           The x axis has the smallest variance, whereas the z axis has the biggest variance.
%   Delta_phi - 3x1 vector with errors in the direction of eigenvectors
%           with respect to the other two directions (Delta_phi is in degrees).
%           Delta_phi(1)==\Delta\phi_12 (between max.var and intermediate var 
%               directions);
%           Delta_phi(2)==\Delta\phi_13 (between max.var and normal
%               directions);
%           Delta_phi(3)==\Delta\phi_23 (between intermediate var and
%               normal directions);
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================
%MVA consists in constructing a matrix M_\mu\nu and solving its eigen-problem.

%% 0. Select data
[selected_times, Bselected] = gap_select_time_range(times, time_range, B);

%% 1. Construct the matrix
% Matrix is M = M1 - M2 = <B_\mu B_\nu> - <B_\mu> <B_\nu>
M1 = (Bselected' * Bselected) / size(Bselected,1);  %size(Bselected,1) === N, number of vectors
mean_B = mean(Bselected,1);
M2 = mean_B' * mean_B;
M = M1 - M2;

%% 2. Solve the eigenvalue problem and determine its min and max variance.
[V, D] = eig(M);
variance = diag(D);
% sort eigenvalues, in case they are not already sorted (as it was in older versions)
[variance,ind] = sort(variance,'ascend');  % store the indices of which 
                                           % columns the sorted eigenvalues come from
RotMatrix = V(:,ind);  % arrange the columns of eigenvectors in this order

% The following lines were taken from MinVariance.h from QSAS:
% // check the positive sense of the eigenvector associated with
% // the minimum eigenvalue and if necessary reverse the signs
% // of all three components i.e. vint x vmax = vector perpendicular
% // to both vint and vmax. If vmin is perpendicular to vint
% // and vmax (i.e. system is ortogonal), then (vinx x vmax).vmin=1
% // If < 0 then not a right-handed system, so reverse signs of vmin

if (dot(cross(RotMatrix(:,2), RotMatrix(:,3)), RotMatrix(:,1)) <0)
    RotMatrix(:,1) = -1 * RotMatrix(:,1);
end

% The following lines were adapted from MinVariance.h from QSAS:
% // ensure that the sign of the average maximum variance component
% // is positive by computing average of sum of B-field data.max
% // eigenvector. If less than zero negate both vmax and vint (so
% // preserving right-handed system)
if (mean(B * RotMatrix(:,3)) < 0)
    RotMatrix(:,2:3) = -1 * RotMatrix(:,2:3);
end

normal = (RotMatrix(:,1))';  % make it also a standard row vector.

%% 3. Compute B in the new frame
if (nargout > 3)
    B_mva = B * RotMatrix;
end

%% 4. Number of points
if (nargout > 4)
    nr_points = length(selected_times);
end

%% 5. Errors using analytical formula
if (nargout > 5)
    % Compute errors using formula (8.23) from chapter 8 of the ISSI'1998
    % book. 
    Delta_phi = zeros(3,1);  % we will use only upper diagonal
    pairs = [1, 2; 1, 3; 2, 3];
    for p=1:length(pairs)
        ci = pairs(p,1);
        cj = pairs(p,2);
        Delta_phi(p) = sqrt((variance(1) / (nr_points-1)) ...
            * (variance(ci) + variance(cj) - variance(1)) / (variance(ci)-variance(cj))^2 );
    end
    Delta_phi = Delta_phi / pi * 180.0;
end

end
