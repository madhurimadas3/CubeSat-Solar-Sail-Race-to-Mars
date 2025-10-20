function  [dot_x] = ODEs(t,x,const,const_moon,flag);
%ODES  
% [dot_x] =ODES(t,x,const) returns x_dot = f(x,t) by specifying the 
% differential equations of the system in first-order form.
%
% INPUT PARAMETERS:
% t = time
% x = system states
% const = a structure that contains relevant physical parameters
% const = a structure that contains relevant physical parameters for the
% Moon
% flag = a flag to determine whether to output x_dot or terms for
% postprocessing
%
% OUTPUT PARAMETERS:
% dot_x = the first-order differential equation evaluated at x and t
%
% Ryan Caverly
% Updated September 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Constants

mu1 = const.mu1;
J2 = const.J2;
Re = const.Re;
kappa = const.beta*const.mu_Sun/const.r_Sun^2;

eye_3 = [0;0;1]; % third column of the identity matrix

% Define the time-varying direction of the Sun
sun_dir = [cos(2*pi/(3600*24*365)*t);sin(2*pi/(3600*24*365)*t);0];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Relevant Positions and Velocities

% First, extract states in a convenient form. 
r_a = x(1:3);
dot_r_a = x(4:6);
r_a_moon = x(7:9);
dot_r_a_moon = x(10:12);
alpha_star_filt = x(13);

% Relative Position to Moon
r_a_rel_moon = r_a - r_a_moon;
r_a_rel_Sun = r_a + sun_dir*const.r_Sun;
r_moon_rel_Sun = r_a_moon + sun_dir*const.r_Sun;

% Norm of position vectors
r_norm = sqrt(r_a'*r_a);
r_norm_moon = sqrt(r_a_moon'*r_a_moon);
r_norm_rel_moon = sqrt(r_a_rel_moon'*r_a_rel_moon);
r_norm_rel_Sun = sqrt(r_a_rel_Sun'*r_a_rel_Sun);
r_norm_moon_rel_Sun = sqrt(r_moon_rel_Sun'*r_moon_rel_Sun);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define McInnes' Optimal Steering Law

% Compute angle and vectors needed for McInnes' steering law
psi_vec = cross(sun_dir,dot_r_a)/norm(cross(sun_dir,dot_r_a));
h_vec = cross(r_a,dot_r_a)/norm(cross(r_a,dot_r_a));

% Correct the signs of psi
if psi_vec'*h_vec > 0
    psi = acos(dot_r_a'*sun_dir/(norm(dot_r_a)*norm(sun_dir)));
else
    psi = -acos(dot_r_a'*sun_dir/(norm(dot_r_a)*norm(sun_dir)));
end

% Compute SIA from McInnes' locally-optimal steering law
alpha_star = 0.5*(psi-asin(sin(psi)/3));

% Place constraints on SIA from McInnes' steeting law
if alpha_star < -const.max_SIA
    alpha_star = -const.max_SIA;
elseif alpha_star > const.max_SIA
    alpha_star = const.max_SIA;
end

% Maximize escape velocity if beyond 1,000,000 km from Earth
if r_norm > 10e8
    alpha_star = 0;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Forces/Accelerations

% Compute solar sail SRP acceleration
Thrust_n = kappa*(cos(alpha_star))^2;

% Compute solar sail normal direction
Normal = (AxisAngle2DCM(psi_vec,abs(alpha_star_filt)))'*sun_dir;

% Compute solar sail acceleration vector
Thrust = Thrust_n*Normal;

% Inverse-Squared Gravitational forces
fg_a = -mu1/r_norm^3*r_a;
fg_a_moon = -mu1/r_norm_moon^3*r_a_moon;

% Additional gravitational force due to J2 perturbation
fg_J2_a = 3*mu1*J2*Re^2/(2*r_norm^5)*((5*r_a(3)^2/r_norm^2-1)*r_a - 2*r_a(3)*eye_3);

% Additional gravitational force due to Moon
fg_moon_pert = const_moon.mu1*(-r_a_rel_moon/r_norm_rel_moon^3 - r_a_moon/r_norm_moon^3);

% Additional gravitational force due to Sun
fg_Sun_pert = const.mu_Sun*(-r_a_rel_Sun/r_norm_rel_Sun^3 + sun_dir/const.r_Sun^2);

% Additional gravitational force on Moon due to Sun
fg_moon_Sun_pert = const.mu_Sun*(-r_moon_rel_Sun/r_norm_moon_rel_Sun^3 + sun_dir/const.r_Sun^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form dot_x = f(x,u) system.

ddot_r_a = fg_a + const.J2_bool*fg_J2_a ...
            + const.moon_bool*fg_moon_pert ...
            + const.Sun_bool*fg_Sun_pert + Thrust;
ddot_r_a_moon = fg_a_moon + const.Sun_bool*fg_moon_Sun_pert;
alpha_filt_dot = -const.max_SIA_rate/const.max_SIA/2*alpha_star_filt + const.max_SIA_rate/const.max_SIA/2*alpha_star;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exit Flag Used for Post Processing

if flag == 0
    dot_x = [dot_r_a;ddot_r_a;dot_r_a_moon;ddot_r_a_moon;alpha_filt_dot];
else
    output.psi = psi;
    output.alpha_star = alpha_star_filt;
    output.alpha_star_rate = alpha_filt_dot;
    output.Thrust = Thrust;
    dot_x = output;
end



