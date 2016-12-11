function [normal, RotMatrix, variance, B_mva, Delta_phi] = gap_mva(B)
%Performs the unconstrained Minimum Variance Analysis (MVA) on the time series B.
%USAGE:
%   [n, RotMatrix, variance, B_mva] = gap_mva(B)
%WHERE
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
%   B_mva - optional array Nx3 with the vectors B transformed into the MVA frame. 
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

%% REFERENCES
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

%MVA consists in constructing a matrix M_\mu\nu and solving its eigen-problem.

%% 1. Construct the matrix
% Matrix is M = M1 - M2 = <B_\mu B_\nu> - <B_\mu> <B_\nu>
M1 = (B' * B) / size(B,1);  %size(B,1) === N, number of vectors
mean_B = mean(B,1);
M2 = mean_B' * mean_B;
M = M1 - M2;

%% 2. Solve the eigenvalue problem and determine its min and max variance.
[V, D] = eig(M);
variance = diag(D);
% sort eigenvalues, in case they are not already sorted (as it was in older versions)
[variance,ind] = sort(variance,'ascend');  % store the indices of which 
                                           % columns the sorted eigenvalues come from
RotMatrix = V(:,ind);  % arrange the columns of eigenvectors in this order
normal = (RotMatrix(:,1))';  % make it also a standard row vector.

%% 3. Compute B in the new frame
if (nargout > 3)
    B_mva = B * RotMatrix;
end

%% 4. Errors using analytical formula
if (nargout > 4)
    % Compute errors using formula (8.23) from chapter 8 of the ISSI'1998
    % book. 
    Delta_phi = zeros(3,1);  % we will use only upper diagonal
    pairs = [1, 2; 1, 3; 2, 3];
    for p=1:length(pairs)
        ci = pairs(p,1);
        cj = pairs(p,2);
        Delta_phi(p) = sqrt((variance(1) / (size(B,1)-1)) ...
            * (variance(ci) + variance(cj) - variance(1)) / (variance(ci)-variance(cj))^2 );
    end
    Delta_phi = Delta_phi / pi * 180.0;
end

end
