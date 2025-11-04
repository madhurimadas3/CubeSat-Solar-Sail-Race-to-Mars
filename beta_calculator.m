close all; clear; clc;

%% Peyton Kramlich - AEM 4331
%% Calculating Solar Sail Lightness Number (Beta)
%% HOW TO USE: change length of boom (l) diagonal of BUS

% input length of boom to determine sail area
l = 7;                    % m, length of boom (INPUTTED VALUE)
diag = sqrt(2*0.2263^2);  % m, diagonal of BUS (INPUTTED VALUE) (currently 12U config)
d = 2*diag;               % m, dist of sail not on boom
D = l-d;                  % m, dist of sail on boom (base
c = sqrt(2*D^2);          % m, pythagoran theorem
A = c^2 - d^2;            % approx area of sail
fprintf("Area of Sail: %0.4f m^2\n", A);

% define space constants
G = 6.6743e-11;         % m^3/kgs
m_sun = 1.989e+33;      % g
m_earth = 5.972e+27;    % g
m_mars = 6.39e+26;      % g

% solar sail materials densities - Goodfellow PEN
rho_sail = 1360;        % kg/m^3
rho_Al = 2700;          % kg/m^3
rho_Cr = 7190;          % kg/m^3
rho_cf = 1600;          % ,, number from google

t_sail = 1.3*10e-6;   % meters
t_Al = 0.1*10e-6;     % meters
t_Cr = 0.015*10e-6;   % meters

% satellite mass calculation (in grams)
% m = 16000; % g      % mass of acs3 cubesat
% m = 14000;          % cubesat bus min mass config
m_camera = 40;
m_panel = 4*289;
m_ais = 375;
m_comms = 52.85;
m_battery = 335;
m_adcs = 1520;
m_avionics = m_camera + m_panel + m_ais + m_comms + m_battery+m_adcs;

m_bus = 800;
m_deploy = 7700; % ACS3 Sail-Boom Subsystem
% m_booms = 4*21.15; % for a 1m long boom (NO)
m_booms = 4*rho_cf*l*(pi*(33e-3 - (33e-3 - 0.12e-3))); % 4 booms, cf density, l, cross-sect area hollow circle
m_structure = m_bus + m_deploy + l*m_booms;

m_sail = A*(rho_sail*t_sail + rho_Al*t_Al + rho_Cr*t_Cr);
m = m_avionics + m_structure + m_sail;
fprintf("Mass of CubeSat: %0.4f g\n", m);

% more constants and equations for beta calculation
rho = 0.88;         % approximate reflectivity coefficient for Al
ab = 1 - rho;       % absorption coefficient

R0 = 1.496e+11;      % m        1 AU, distance from Sun to Earth
P0 = 4.563;          % N/m^2    SRP at Earth
Rm = 2.28e+11;       % m        distance from Sun to Mars
Pm = P0 * (R0/Rm)^2; % N/m^2    SRP at Mars
R = 1;               % m        distance from Sun to sail
P = P0 * (R0/R)^2;   % N/m^2    local SRP

% calculate area-to-mass ratios (sigma)
sigma_sail = m/A;
sigma_earth = (2*P0*R0^2)/G/m_sun;
sigma_mars = (2*Pm*Rm^2)/G/m_sun;

% calculate beta
beta = sigma_earth/sigma_sail;
fprintf("Beta: %.4f \n", beta);
