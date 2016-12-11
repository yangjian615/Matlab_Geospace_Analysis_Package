function [remainder] = gap_remainder(A, P)
%Computes remainder from a division.
%USAGE:
%   [remainder] = gap_remainder(A, P)
%WHERE
%   1) A is a matrix and P is a scalar -- computes remainder of all elements of A by the same P.
%   2) A and P are matrices of the same size -- computes remainder elementweise. 
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

remainder = A - (fix(A./P).*P);

end
