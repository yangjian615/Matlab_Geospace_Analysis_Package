function C = gap_cross(A, B)
%The cross product of vectors, provided that they are in the row form.
%   The Matlab's included cross() function does too many checks and rotations, which slows
%   down greatly the computation. We can skip those checks by trusting that scripts that
%   call this function provide arrays with correct dimensionality and are row vectors. 
%USAGE:
%   C = gap_cross(A, B)
%WHERE
%   A, B - are Nx3 array with data of cartesian type
%OUTPUT
%   C   - array Nx3 where C(i,:) = A(i,:) x B(i,:)
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

C = [A(:,2) .* B(:,3) - A(:,3) .* B(:,2), ...
     A(:,3) .* B(:,1) - A(:,1) .* B(:,3), ...
     A(:,1) .* B(:,2) - A(:,2) .* B(:,1)];

end
