function [m, com, I] = Mass_Properties_scaled(s)
%MASS_PROPERTIES_SCALED  Link masses/inertias for a uniformly scaled robot.
%
%   [m, com, I] = Mass_Properties_scaled(s)
%
% Same body layout and geometry as the nominal Mass_Properties() local
% function in rabbit_energy_model_generalized_Lagrange.m (torso, two thighs,
% two shins), but every link mass and every entry of its 4x4 pseudo-inertia
% matrix is multiplied by s. That is the Lagrangian-derivation counterpart of
% Chapter 4's mass_scale perturbation: scaling every link's density by s
% leaves its shape and COM unchanged and scales both the COM-frame rotational
% inertia and the parallel-axis correction m*(...) by the same factor s,
% since both are linear in mass for a fixed geometry. Passing s = 1 must
% reproduce Mass_Properties() exactly.
%
% Input
%   s : mass_scale, positive finite scalar
%
% Outputs
%   m   : 5 x 1   [torso thigh1 shin1 thigh2 shin2] masses, scaled by s
%   com : 4 x 5    homogeneous COM positions (unchanged by s)
%   I   : 4 x 4 x 5   pseudo-inertia matrices, scaled by s
%
% See also RABBIT_ENERGY_MODEL_GENERALIZED_LAGRANGE, RABBIT_GENERATE_CASE_DYNAMICS.

if ~(isscalar(s) && isfinite(s) && s > 0)
    error('Mass_Properties_scaled:s', ...
          'mass_scale must be a positive finite scalar (got %s).', mat2str(s));
end

% +++++++ torso +++++++
Mt = 47; %Kg
X = 0; Y = -0.75/2; Z = 0;
comt = [X Y Z]';
Ixx = 0.94; Ixy = 0; Ixz = 0;
Iyx = 0;    Iyy = 0.235; Iyz = 0;
Izx = 0;    Izy = 0;     Izz = 0.94;
Ixx = Ixx + Mt*(Y^2+Z^2);
Iyy = Iyy + Mt*(X^2+Z^2);
Izz = Izz + Mt*(Y^2+X^2);
I_t = [ (-Ixx+Iyy+Izz)/2   Ixy               Ixz               Mt*X; ...
        Iyx                (Ixx-Iyy+Izz)/2   Iyz               Mt*Y; ...
        Izx                Izy               (Ixx+Iyy-Izz)/2   Mt*Z; ...
        Mt*X               Mt*Y              Mt*Z              Mt];

% +++++++ thigh 1 +++++++
M1 = 10; %Kg
X = 0; Y = 0.5/2; Z = 0;
com1 = [X Y Z]';
Ixx = 0.2; Ixy = 0; Ixz = 0;
Iyx = 0;   Iyy = 0.04; Iyz = 0;
Izx = 0;   Izy = 0;    Izz = 0.2;
Ixx = Ixx + M1*(Y^2+Z^2);
Iyy = Iyy + M1*(X^2+Z^2);
Izz = Izz + M1*(Y^2+X^2);
I_1 = [ (-Ixx+Iyy+Izz)/2   Ixy               Ixz               M1*X; ...
        Iyx                (Ixx-Iyy+Izz)/2   Iyz               M1*Y; ...
        Izx                Izy               (Ixx+Iyy-Izz)/2   M1*Z; ...
        M1*X               M1*Y              M1*Z              M1];

% +++++++ shin 1 +++++++
M2 = 3.5; %Kg
X = 0; Y = 0.5/2; Z = 0;
com2 = [X Y Z]';
Ixx = 0.07; Ixy = 0; Ixz = 0;
Iyx = 0;    Iyy = 0.014; Iyz = 0;
Izx = 0;    Izy = 0;     Izz = 0.07;
Ixx = Ixx + M2*(Y^2+Z^2);
Iyy = Iyy + M2*(X^2+Z^2);
Izz = Izz + M2*(Y^2+X^2);
I_2 = [ (-Ixx+Iyy+Izz)/2   Ixy               Ixz               M2*X; ...
        Iyx                (Ixx-Iyy+Izz)/2   Iyz               M2*Y; ...
        Izx                Izy               (Ixx+Iyy-Izz)/2   M2*Z; ...
        M2*X               M2*Y              M2*Z              M2];

% +++++++ thigh 2 +++++++
M3 = 10; %Kg
X = 0; Y = 0.5/2; Z = 0;
com3 = [X Y Z]';
Ixx = 0.2; Ixy = 0; Ixz = 0;
Iyx = 0;   Iyy = 0.04; Iyz = 0;
Izx = 0;   Izy = 0;    Izz = 0.2;
Ixx = Ixx + M3*(Y^2+Z^2);
Iyy = Iyy + M3*(X^2+Z^2);
Izz = Izz + M3*(Y^2+X^2);
I_3 = [ (-Ixx+Iyy+Izz)/2   Ixy               Ixz               M3*X; ...
        Iyx                (Ixx-Iyy+Izz)/2   Iyz               M3*Y; ...
        Izx                Izy               (Ixx+Iyy-Izz)/2   M3*Z; ...
        M3*X               M3*Y              M3*Z              M3];

% +++++++ shin 2 +++++++
M4 = 3.5; %Kg
X = 0; Y = 0.5/2; Z = 0;
com4 = [X Y Z]';
Ixx = 0.07; Ixy = 0; Ixz = 0;
Iyx = 0;    Iyy = 0.014; Iyz = 0;
Izx = 0;    Izy = 0;     Izz = 0.07;
Ixx = Ixx + M4*(Y^2+Z^2);
Iyy = Iyy + M4*(X^2+Z^2);
Izz = Izz + M4*(Y^2+X^2);
I_4 = [ (-Ixx+Iyy+Izz)/2   Ixy               Ixz               M4*X; ...
        Iyx                (Ixx-Iyy+Izz)/2   Iyz               M4*Y; ...
        Izx                Izy               (Ixx+Iyy-Izz)/2   M4*Z; ...
        M4*X               M4*Y              M4*Z              M4];

m   = s * [Mt M1 M2 M3 M4]';
com = [comt com1 com2 com3 com4];
com = [com; ones(1,5)];

I(:,:,1) = s * I_t;
I(:,:,2) = s * I_1;
I(:,:,3) = s * I_2;
I(:,:,4) = s * I_3;
I(:,:,5) = s * I_4;

end
