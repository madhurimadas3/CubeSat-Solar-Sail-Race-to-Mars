% Constants file.

% Define Solar Sail Parameters
const.beta = 0.025; % ACS3 (12U) = 0.008, HiPERSail (27U) = 0.025
const.max_SIA = 65*pi/180; % (rad)
const.max_SIA_rate = 0.02*pi/180; % (rad/s) (Solar Cruiser used 0.02 deg/s max rate)

% adding this part (MD)
% Additional spacecraft properties
const.mass = 12;        % kg (example for 12U sail)
const.A = 100;           % m^2 (sail area) -- change this as necessary
const.Cr = 1.8;         % reflectivity coefficient
const.c = 3e8;          % m/s (speed of light)
const.S0 = 1367;        % W/m^2 (solar constant)
const.AU = 1.496e11;    % m (1 AU)


% Perturbation Options (1: on, 0: off)
const.J2_bool = 1; % Effect of J2 perturbation
const.moon_bool = 1; % Effect of Moon gravitational pull on solar sail
const.Sun_bool = 1; % Effect of Sun gravitational pull on solar sail

% Gravitational constants
const.mu1 = 3.986e14; % Earth gravitational constant (m^3/s^2)
const_moon.mu1 = 0.049e14; % Earth gravitational constant (m^3/s^2)
const.J2 = 1.08262645e-3; % J2 perturbation coefficient (set to 0 if no J2 perturbation)

% Other constants
const.Re = 6371.2e3; % radius of the Earth (m)
const.mu_Sun = 1.327e20;
const.r_Sun = 149e9; %(m)

% Define Initial Solar Sail Orbit
const.a = 42464e3; % semimajor axis of orbit (m)
const.e = 0.0002771; % eccentricity of orbit
const.inc = 23.5*pi/180; % inclination of orbit (rad)
const.omega = rand*180*pi/180; % argument of perigeee (rad)
const.Omega = rand*180*pi/180; % right ascension of the ascending node (rad)
const.t0 = 0; % time of perigee passage (s)

% Define Moon's orbit
const_moon.a = 383200e3; % semimajor axis of orbit (m)
const_moon.e = 0.063; % eccentricity of orbit
const_moon.inc = 5.1*pi/180; % inclination of orbit (rad)
const_moon.omega = rand*180*pi/180; % argument of perigeee (rad)
const_moon.Omega = rand*180*pi/180; % right ascension of the ascending node (rad)
const_moon.t0 = 0; % time of perigee passage (s)




