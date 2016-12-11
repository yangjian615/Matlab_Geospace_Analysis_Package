function [curlB, K, divB, GradB] = gap_curlometer(pos_vec4s, B_vec4s)
%Computes the curl of the magnetic field provided by four satellites. 
%USAGE:
%   [curlB, K, divB, GradB] = gap_curlometer(pos_vec4s, B_vec4s)
%WHERE
%   pos_vec4s - is a cell-array with 4 cells, corresponding to position data from the four
%               satellites.
%   B_vec4s - is a cell-array with 4 cells, corresponding to magnetic field data from the four
%               satellites. 
%OUTPUT
%   curlB - the estimated curl of the magnetic field. The unit of this estimation is equal
%               to units of B_vec4s divided by units of pos_vec4s.
%               Cluster's FGM data have position in units of 'km' and the magnetic field in
%               units of 'nT'. Thus, the units of curlB are 'nT/km'. The following expression
%               can be used to convert curlB into current densities of units A/m^2:
%                   J = curlB * 1e-12 / (1.25663706e-06);
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

Nsat = length(B_vec4s);
if (Nsat == 4)
    
    %% 0. Check if the geometry is good enough.
    
    %% 1. Compute reciprocal vectors
    Nmeas = size(B_vec4s{1},1);
    K = zeros(Nmeas, 3, Nsat);
    B_adapted4s = zeros(Nmeas, 3, Nsat);
    c = [1, 2, 3, 4, 1, 2, 3, 4];  % cyclic permutation for satellites
    for sat = 1:Nsat
        S = gap_cross( (pos_vec4s{c(sat+2)}(:,:) - pos_vec4s{c(sat+1)}(:,:)) ...
                       , (pos_vec4s{c(sat+3)}(:,:) - pos_vec4s{c(sat+1)}(:,:)) );
        Vol = dot( (pos_vec4s{c(sat+0)}(:,:) - pos_vec4s{c(sat+1)}(:,:)), S , 2);
%         Volume_minim = vnorm(S, 2) .* ...
%             vnorm((pos_vec4s{c(sat+0)} - pos_vec4s{c(sat+1)}), 2)/ 5/3;
%         Volume_minim = ones(length(Volume_minim), 1) * 1.0e6;  
        %uncomment above lines in order to ignore the effect of min volume
        
        Volume_minim = 0.001;
        index = find(abs(Vol) >= Volume_minim);
        for x = 1:3
            K(index,x,sat) = (S(index,x) ./ Vol(index));
        end
        B_adapted4s(:,:,sat) = B_vec4s{sat};
    end
    clear sat c x index
    
    %% 2. Compute the gradient
    GradB = zeros(3, 3, Nmeas);  %the gradient tensor of a vector
    for t = 1:Nmeas
        GradB(:,:,t) = squeeze(K(t,:,:)) * ((squeeze(B_adapted4s(t,:,:)))');
    end
    
    %% 3. Compute curl B
    curlB = zeros(Nmeas, 3);
    c = [1, 2, 3, 1, 2, 3];  % cyclic permutation for vector components
    for x = 1:3
        curlB(:,x) = GradB(c(x+1), c(x+2),:) - GradB(c(x+2), c(x+1),:);
    end
    
    %% 4. Compute divergence
    if (nargout > 2)
        divB = zeros(Nmeas,1);
        for t=1:Nmeas
            divB(t) = trace(squeeze(GradB(:,:,t)));
        end
    end
else
    error('MATLAB:gap_curlometer:InvalidNrSats',...
        'Currents are not computed. Curlometer requires magnetic field from 4 satellites.');
end

end
