function [angle] = gap_angle_between_vectors(vect1, vect2)
%Computes angle between vectors
%Usage:
%   [angle] = gap_angle_between_vectors(vect1, vect2)
%Where
%   vect1, vect2 - array Nx3 of N vectors.
%   angle - array Nx1 of angles in the range 0..180 degrees.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

angle = real(acosd(dot(vect1, vect2, 2) ./ vnorm(vect1,2) ./ vnorm(vect2,2)));

end
