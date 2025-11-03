% Created by Ryan James Caverly
% September 2025

% Modified by Madhurima Das
% Date: Oct 29, 2025
% Simulation of a solar sail deployed from a GEO graveyard orbit performing
% an orbit-raising maneuver modified from McInnes' locally optimal steering
% law

clear all
clc
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants 
constants_struct % All constants in one file. 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (MD)
% Lauch Date
ref_launch_date = datetime(2026,1,1,0,0,0); % reference epoch (t=0) % CHANGE THIS TO ACTUAL LAUNCHDATE
launch_date_offset_days = 0; % CHANGE THIS AS NEEDED/for alignment testing
launch_date = ref_launch_date + days(launch_date_offset_days);
launch_date_offset_sec = seconds(launch_date - ref_launch_date);

fprintf('Simulated launch date: %s\n', datestr(launch_date));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions.

% Find initial position and velocity given orbital elements
[r_a0,v_a0] = findRandV(const.mu1,const.e,const.inc,const.a,const.omega,const.Omega,const.t0,launch_date_offset_sec);
[r_a0_moon,v_a0_moon] = findRandV(const.mu1,const_moon.e,const_moon.inc,const_moon.a,const_moon.omega,const_moon.Omega,const_moon.t0,launch_date_offset_sec);

% Computations to Properly Initialize alpha_star for SIA rate filter
sun_dir = [1;0;0];
psi_vec = cross(sun_dir,v_a0)/norm(cross(sun_dir,v_a0));
h_vec = cross(r_a0,v_a0)/norm(cross(r_a0,v_a0));
if psi_vec'*h_vec > 0
    psi = acos(v_a0'*sun_dir/(norm(v_a0)*norm(sun_dir)));
else
    psi = -acos(v_a0'*sun_dir/(norm(v_a0)*norm(sun_dir)));
end
alpha_star0 = 0.5*(psi-asin(sin(psi)/3));
if alpha_star0 < -const.max_SIA
    alpha_star0 = -const.max_SIA;
elseif alpha_star0 > const.max_SIA
    alpha_star0 = const.max_SIA;
end

IC = [r_a0;v_a0;r_a0_moon;v_a0_moon;alpha_star0];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time.

sim_length_days = 1*365; % simulation length in days % 3years is a good time to escape

t0 = 0; % s
t_max = 3600*24*sim_length_days; % s
t_div = 20001;  % number of steps to divide the time series into.
t_span = linspace(t0,t_max,t_div); % Total simulation time.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation options.
 options = odeset('AbsTol',1e-10,'RelTol',1e-10); % This changes the integration tolerence.

flag = 0;
[t,x_out] = ode78(@ODEs,t_span,IC,options,const,const_moon,flag);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing
post_processing

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data
plot_script_v1
