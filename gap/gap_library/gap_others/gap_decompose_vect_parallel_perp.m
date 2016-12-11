function [paral_abs, paral_vec, perp_vec] = gap_decompose_vect_parallel_perp(in_vect, ref_vect)
%Decomposes vector into parallel and perpendicular components with respect 
%   to a reference vector (e.g. background magnetic field).
%Usage:
%   [paral_abs, paral_vect, perp] = gap_decompose_vect_parallel_perp(in_vect, ref_vect)
%Where:
%   in_vect and ref_vect are Nx3 matrices
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

unit_ref_vect = ref_vect ./ repmat(vnorm(ref_vect, 2), 1, 3);
paral_abs = dot(in_vect, unit_ref_vect, 2);
if (nargout > 1)  %show also the parallel component as a full vector
    paral_vec = repmat(paral_abs, 1,3) .* unit_ref_vect;
    if (nargout > 2)  %show also the perpendicular component
        perp_vec = in_vect - paral_vec;
    end
end

end
