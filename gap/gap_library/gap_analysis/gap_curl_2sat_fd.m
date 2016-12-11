function [times_Jfd, Jfd, u, v, z, r4, p4, L, l, M, m, eps, lambda, mu, MM] ...
    = gap_curl_2sat_fd(j_times ...
    , B2joined, pos2joined, assumption, estim_e, jump_displace)
%Estimates the curl of the magnetic field, based on data from two
%   satellites that are flying quasi-parallel. 
%USAGE:
%   [times_Jfd, Jfd, u, v, z, r4, p4, L, l, M, m, eps, lambda, mu, MM] ...
%    = gap_curl_2sat_fd(j_times, B2joined, pos2joined, assumption, estim_e, jump_displace)
%WHERE
%   j_times    - Nx1 array with time tags, where N is the nr of
%               measurements and each time tag is in the datenum format
%               (i.e. number of days since 2000-Jan-01). 
%   B2joined   - is a cell-array with 2 cells, corresponding to magnetic 
%               field data from the two satellites. Each cell is Nx3
%               array. The data has to be already time joined before
%               calling this function.
%   pos2joined - is a cell-array with 2 cells, corresponding to position 
%               data from the two satellites. Each cell is Nx3 array.
%   assumption - is one of options to estimate the component of the
%               gradient, which is perpendicular to the plane formed by the
%               two satellites and their velocity vector. 
%               'no assum'  - means that the normal component of the
%                   gradient will not be estimated. The result will contain
%                   only the gradient in the plane of the two satellites.
%               'perp'      - gradient is assumed to be perpendicular to a
%                   given vector 'estim_e'.
%               'paral'     - gradient is parallel to the vector 'estim_e'.
%               'force free'- assumes that the curl of the field is
%                   parallel to the background field (j x B = 0). The
%                   background field is given in the variable 'estim_e'.
%   estim_e - array 1x3 representing the vector used in various
%               assumptions. See the 'assumption' variable.
%   jump_displace - is an indicator for the time interval:
%               \Delta t (seconds) = (j_times(1+jump_displace) -
%                                     j_times(1)) * 24.0 * 3600.0
%               It is assumed that the gradient is constant at the scales
%               of 2*\Delta t. If it is smaller, then the round-off
%               errors are larger. 
%OUTPUT
%   times_Jfd - Nx1 of datenums. Time tags of the resulting current estimations. 
%   Jfd - Nx3 the estimated curl of the magnetic field. The unit of this 
%               estimation is equal to units of B2joined divided by units 
%               of pos2joined. 
%               Cluster's FGM data have position in units of 'km' and the magnetic field in
%               units of 'nT'. Thus, the units of curlB are 'nT/km'. The following expression
%               can be used to convert curlB into current densities of units A/m^2:
%                   J = curlB * 1e-12 / (1.25663706e-06);
%   u, v, z - unit vectors in the reference system associated with the pair
%               of satellites. v is in the direction of center of mass'
%               velocity. u is perpendicular to v, but in the plane formed
%               by v and the separation of the two satellites. z completes
%               the triad.
%               If the rotation matrix is rotMatr = [u; v; z]', then the
%               transformation between uvz and initial frames is:
%               p_gsm = p_uvz * rotMatr';
%   r4 - coordinates of the measurement points in the (u,v,z) frame:
%               [r_a^+; r_a^-; r_b^+; r_b^-]. See expressions (26) and (27)
%               in the Vogt et al. (2013) paper. 
%   p4 - 4x3 with [pa_plus; pa_minus; pb_plus; pb_minus] in the (u, v, z)
%               frame.
%   L, l, M, m - discretization parameters. See relations (22) to (25) in
%               the paper.
%   eps, lambda, mu - Another way to express above parameters.
%   MM - 3x3xN is a gradient matrix for each of the output time tags.
%
%This program comes with ABSOLUTELY NO WARRANTY. See "Geospace Analysis 
%   Package" project's LICENSE file for more information.

%% REFERENCES
%   Vogt, J., Sorbalo, E., He, M., and Blagau, A. (2013). Gradient 
%   estimation using configurations of two or three spacecraft. Annales 
%   Geophysicae (09927689), 31(11), 1913-1927. doi:10.5194/angeo-31-1913-2013

%% Copyright (c) 2014 "Geospace Analysis Package" project
%   This function is a part of the "Geospace Analysis Package" project.
%   Please see the GAP project's license before modifying or redistributing
%   this file.

%% ========================================================================

if ~iscell(B2joined) || ~iscell(pos2joined)
    error('MATLAB:gap_curl_2sat_fd:wrong_input', 'Input must be cells.');
end

ssat1 = 1; ssat2 = 2;

%% Special units: L, l, M, m
r12 = mean(pos2joined{ssat2} - pos2joined{ssat1}, 1) / 2;
dr_1=(pos2joined{ssat1}((1+jump_displace):end,:) - pos2joined{ssat1}(1:(end-jump_displace),:));
dr_2=(pos2joined{ssat2}((1+jump_displace):end,:) - pos2joined{ssat2}(1:(end-jump_displace),:));

Va = mean(dr_1,1) ./ (repmat(j_times(1+jump_displace)-j_times(1),1,3)*24.0*3600);
Vb = mean(dr_2,1) ./ (repmat(j_times(1+jump_displace)-j_times(1),1,3)*24.0*3600);

V_star = (Va + Vb) / 2;
V_Delta = (Vb - Va) /2;
vv = V_star ./ repmat(vnorm(V_star,2),1,3);
v = vv(1,:);  %take only one. Fix this later

ll = dot(r12,repmat(vv,size(r12,1),1),2);  %multiple
ll = ll(1);
uu = r12(1,:) - repmat(ll,1,3) .* vv;
uu = uu ./ repmat(vnorm(uu,2),1,3);
u = uu(1,:);  %take only one. FIX IT!

%Perpendicular to 'uu' and 'vv':
z = gap_cross(uu, vv);
z = z / vnorm(z,2);

dtime_vogt = mean(j_times(1+jump_displace:end) - j_times(1:(end-jump_displace)))*24*3600;

l = (ll);
L = (dot(r12,repmat(uu,size(r12,1),1),2));  %multiple
L = L(1);
M = (dot(V_star,repmat(vv,size(V_star,1),1),2)) * dtime_vogt;  %one
m = (dot(V_Delta,repmat(uu,size(V_Delta,1),1),2)) * dtime_vogt;  %one

eps = m ./ M;
lambda = l ./ L;
mu = M ./ L;

%Compute p_a^\pm, p_b^\pm and finally gradient
pa_plus = [-1.0/4/L - l/(4*L*M),...
    +1.0 / (4*M), 0];
pa_minu = [-1.0/4/L + l/(4*L*M),...
    -1.0 / (4*M), 0];
pb_plus = [+1.0/4/L - l/(4*L*M),...
    +1.0 / (4*M), 0];
pb_minu = [+1.0/4/L + l/(4*L*M),...
    -1.0 / (4*M), 0];

p4 = [pa_plus; pa_minu ...
    ; pb_plus; pb_minu];

% error_in_B = 1;  %nT
% err_grad = sqrt(sum(vnorm(p,2).^2)) * error_in_B;


%% Compute amplification. Formulas (26) and (27) from the paper
r4 = [  -L - m, -l + M, 0 ; ...
    -L + m, -l - M, 0 ; ...
    L + m,  l + M, 0 ; ...
    L - m,  l - M, 0];

% amplification = trace(p * p') * 0.25 * trace(r4*r4');


%%
N = length(B2joined{ssat1});
Jfd = zeros(N,3);

%Next matrix is to transform into satellite frame
rotMatr = [u; v; z]';
p_gsm = p4 * rotMatr';
B_rot = B2joined;

e_z = [0,0,1] * rotMatr';  %z in the (u,v,z) frame
e = estim_e;
if strcmpi(assumption, 'no assum')
    assumpt_i = 0;
elseif strcmpi(assumption, 'perp')
    % Perpendicular assumption
    assumpt_i = 1;
    e = e / vnorm(e,2);
    dot_e_z = dot(e, e_z, 2);
%     fprintf('Quality angle: %.1f (best if ~0 deg). ' ...
%         , gap_angle_between_vectors(e, e_z));
elseif strcmpi(assumption, 'paral')
    % Parallel assumption
    assumpt_i = 2;
    e = e ./ repmat(vnorm(e,2), 1, 3);
    cross_e_z = vnorm(gap_cross(e, repmat(e_z, size(e,1))),2) .^2;
    dot_e_z = dot(e, repmat(e_z, size(e,1)), 2);
%     fprintf('Quality angle: %.1f (best if ~90 deg). ' ...
%         , gap_angle_between_vectors(e, e_z));
elseif strcmpi(assumption, 'force free')
    assumpt_i = 3;
else
    % assumpt_i = 4;  % default value
    error('MATLAB:gap_curl_2sat_fd:assumptions', 'Assumption not recognized');
end
%--
MM = zeros(3, 3, N);
for t=(1+jump_displace):(N-jump_displace)
    B_t = [
        B_rot{ssat1}(t+jump_displace,:); ...
        B_rot{ssat1}(t-jump_displace,:); ...
        B_rot{ssat2}(t+jump_displace,:); ...
        B_rot{ssat2}(t-jump_displace,:) ];
    
    MM(:,:,t) = p_gsm' * B_t;
    MMz = zeros(1,3);
    %% Estimate normal components:
    if (assumpt_i == 1)
        % Perpendicular assumption
        MMz = - (e * (MM(:,:,t))) / dot_e_z;
    elseif (assumpt_i == 2)
        % Parallel assumption
        MMz =  (e * (MM(:,:,t))) * dot_e_z / cross_e_z;
    elseif (assumpt_i == 3)
        % Force free assumption
        div_p_B = 0;
        cross_p_B = zeros(1,3);
        for point = 1:4
            div_p_B = div_p_B + p4(point,:) * B_t(point,:)';
            cross_p_B = cross_p_B + gap_cross(p4(point,:), B_t(point,:));
        end
        tmp1 = gap_cross(e_z, gap_cross(cross_p_B, e));
        MMz = gap_cross(tmp1, e_z) / dot(e_z, e, 2) - div_p_B * e_z;
    end
    %--
    MM(:,:,t) = MM(:,:,t) + rotMatr(:,3) * MMz;
    
    Jfd(t,1) = MM(2,3,t) - MM(3,2,t);
    Jfd(t,2) = MM(3,1,t) - MM(1,3,t);
    Jfd(t,3) = MM(1,2,t) - MM(2,1,t);
end
t_indices = (1+jump_displace):(size(Jfd,1)-jump_displace);
MM = MM(:,:,t_indices);
Jfd = Jfd(t_indices,:);
times_Jfd = j_times(1+jump_displace:end-jump_displace);

end

