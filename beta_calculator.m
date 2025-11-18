close all; clear; clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peyton Kramlich - AEM 4331
% Created:      September 2025
% Last updated: November 2025
% Calculating Solar Sail Lightness Number (Beta)
% HOW TO USE: change length of boom (l) diagonal of BUS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- initialize variables and calculate beta -----
% input length of boom to determine sail area
L = 7;                 % m, length of boom (INPUTTED VALUE)
% currently 12U config
bus_w = 0.2;           % m, width of the bus (INPUTTED VALUE)
bus_L = bus_w;            % m, length of the bus (INPUTTED VALUE)
bus_h = 0.400263;         % m, height of the bus (INPUTTED VALUE)
diag = sqrt(bus_L^2 + bus_w^2); % m, diagonal of BUS
d = 2*diag;               % m, dist of sail not on boom
D = L-d;                  % m, dist of sail on boom (base
C = sqrt(2*D^2);          % m, pythagorean theorem for long side of sail
c = sqrt(2*d^2);          % m, pythagorean theorem for short side of sail
h = sqrt(D^2-((C-c)/2)^2);% m, height of sail
A = 4*0.5*(C+c)*h;        % m^2, approx total area of sail

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
rho_cf = 1800;          % kg/m^3, number from google

t_sail = 1.3e-6;   % meters
t_Al = 0.1e-6;     % meters
t_Cr = 0.015e-6;   % meters

% satellite mass calculation (in grams)
m_camera = 40;
m_panel = 4*289;
m_ais = 375;
m_comms = 52.85;
m_battery = 670;
m_adcs = 1520;
m_avionics = m_camera + m_panel + m_ais + m_comms + m_battery + m_adcs;

m_bus = 772;
m_cable_motors = 4*79.8;
m_deploy_motors = 40;
m_sail_spool = 4*211;
m_deploy_mech = 7700;
%m_deploy = 7700; % ACS3 Sail-Boom Subsystem
m_deploy = m_cable_motors+m_deploy_motors+m_sail_spool+m_deploy_mech;

% 4 booms, cf density, l, cross-sect area hollow circle
od = 33e-3;       % m, outer diameter
t = 0.12e-3;      % m, wall thickness
or = od/2;        % m, outer radius
ir = or - t;      % m, inner radius
A_cf = pi*(or^2 - ir^2);    % m^2, area of carbon fiber
m_booms = 4*rho_cf*L*A_cf*1000;
m_structure = m_bus + m_deploy + m_booms;

m_sail = A*(rho_sail*t_sail + rho_Al*t_Al + rho_Cr*t_Cr)*1000;
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- plot the sail -----

% boom length along x and y-axes
x_boom = -L:0.01:L;
y_boom = -L:0.01:L;

% square bus centered at origin
x_bus = [bus_w/2  bus_w/2 -bus_w/2 -bus_w/2  bus_w/2];
y_bus = [bus_w/2 -bus_w/2 -bus_w/2  bus_w/2  bus_w/2];

% solar panels deployed attached to side of bus
% top
x_panel1 = [bus_w/2  bus_w/2       -bus_w/2       -bus_w/2  bus_w/2];
y_panel1 = [bus_w/2  bus_w/2+bus_h  bus_w/2+bus_h  bus_w/2  bus_w/2];

% left
x_panel2 = [-bus_w/2  -bus_w/2-bus_h  -bus_w/2-bus_h  -bus_w/2  -bus_w/2];
y_panel2 = [ bus_w/2    bus_w/2       -bus_w/2        -bus_w/2   bus_w/2];

% bottom
x_panel3 =    [bus_w/2  bus_w/2       -bus_w/2       -bus_w/2  bus_w/2];
y_panel3 = -1*[bus_w/2  bus_w/2+bus_h  bus_w/2+bus_h  bus_w/2  bus_w/2];

% right
x_panel4 = -1*[-bus_w/2  -bus_w/2-bus_h  -bus_w/2-bus_h  -bus_w/2  -bus_w/2];
y_panel4 =    [bus_w/2    bus_w/2        -bus_w/2        -bus_w/2   bus_w/2];

% solar sails
x_sail1 = [L  bus_w/2  bus_w/2  d  L]; % northeast quad
y_sail1 = [bus_w/2  L  d  bus_w/2  bus_w/2];

x_sail2 = -1*[L  bus_w/2  bus_w/2  d  L]; % northwest quad
y_sail2 = [bus_w/2  L  d  bus_w/2  bus_w/2];

x_sail3 = -1*[L  bus_w/2  bus_w/2  d  L]; % southwest quad
y_sail3 = -1*[bus_w/2  L  d  bus_w/2  bus_w/2];

x_sail4 = [L  bus_w/2  bus_w/2  d  L]; % southeast quad
y_sail4 = -1*[bus_w/2  L  d  bus_w/2  bus_w/2];

% plot the satellite and sail
figure(1)
plot(x_boom, zeros(numel(x_boom),1), 'Color', '#808080', "LineWidth",2); hold on;
plot(zeros(numel(y_boom),1), y_boom, 'Color', '#808080', "LineWidth",2);
fill(x_bus, y_bus, 'k'); % filled black shape

fill(x_panel1, y_panel1, [211 211 211]/255);
fill(x_panel2, y_panel2, [211 211 211]/255);
fill(x_panel3, y_panel3, [211 211 211]/255);
fill(x_panel4, y_panel4, [211 211 211]/255);

fill(x_sail1, y_sail1, [189 213 231]/255);
fill(x_sail2, y_sail2, [189 213 231]/255);
fill(x_sail3, y_sail3, [189 213 231]/255);
fill(x_sail4, y_sail4, [189 213 231]/255);

ylim([-L-0.1 L+0.5]); axis equal;
xlabel('X (m)'); ylabel('Y (m)');
lgd = legend('Booms', '', 'CubeSat Bus', 'Solar Panels', '', '', '', ...
             'Sail', '', '', '', 'Location', 'northwest');
legendPos = lgd.Position;

str = sprintf('\\beta=%.4f\n L_{boom}=%.2f', beta, L);
annotation('textbox', [legendPos(1), legendPos(2)-0.08, legendPos(3), 0.05], ...
    'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', ...
    'k');

% exportgraphics(figure(1), 'solarSail.png', 'ContentType', 'vector');
