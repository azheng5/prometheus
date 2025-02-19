clear;clc;close all

k = 1.9027e-4;
H = 11100; % m
h0 = 100000;
mu_mars = 4.282837e13; % grav param of mars
R_mars = 3390e3; % radius of mars
r_orbiter = 370000 + R_mars;
rho0 = 0.02; % kg/m^3
g = 3.71; % m/s^2
gamma = [1:0.1:10] .* (pi/180); % flight path angle (rad)
gamma0 = 10 .* (pi/180);
m = 200; % mass (kg)
r_b = 2; % radius of base (m)
A = pi*r_b; % reference area
v_atm = 3300; % reentry velocity (m/s)
r_n = 0.1; % radius of nosecone (m)
C_p_max = 2;
delta_c = 60 * (pi/180); % degrees - assume will change with more design
CD = C_p_max * (0.5 .* (1-sin(delta_c).^4)*(r_n/r_b)^2 + (sin(delta_c)^2)*(1-(r_n/r_b)^2*cos(delta_c)));
beta = m/(CD*A); % ballistic coefficient
C = -(rho0.*H)./(2.*beta.*sin(gamma));

t_span = 0:0.1:5000;
init_state = [h0;v_atm;gamma0;0];
config.beta = beta;
config.rho0 = rho0;
config.H = H;
config.R_mars = R_mars;
config.g = g;
config.k = k;
config.r_n = r_n;
options = odeset('Events',@edl_events);
[time, state] = ode45(@(t,x)edl_eom(t,x,config), t_span, init_state, options);

[h_n_max, n_max, q_dot_s_max, J_s] = compute_performance(gamma0, v_atm);

a = zeros(length(time),1);
q_dot = zeros(length(time),1);
for i = 1:length(time)
    [xdot,log] = edl_eom(time(i),state(i,:)',config);
    a(i) = log.a;
    q_dot(i) = log.q_dot;
end

figure(1)
plot(time,state(:,1))
xlabel('Time (sec)')
ylabel('Altitude (m)')

figure(2)
plot(time,state(:,2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')

figure(3)
plot(time,state(:,3) .* (180/pi))
xlabel('Time (sec)')
ylabel('Flight path angle (deg)')

figure(4)
plot(time,a)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')

figure(5)
plot(time,q_dot)
xlabel('Time (sec)')
ylabel('Heat rate')

%% Compute thickness
% v_atm = [1600:1:3600]; % reentry velocity (m/s)
[h_n_max, n_max, q_dot_s_max, J_s] = compute_performance(10*(pi/180), 3300);
[q_dot_s_max_sim index] = max(q_dot);
q_dot_s_max_alt = state(index,1);
q_dot_s_max_vel = state(index,2); % m/s

emissivity = 0.9; % PICA - 0.9; SLA - 0.7
stefanBoltzmannConst = 5.67e-8;

if q_dot_s_max_alt < 7000
    T = -31  - 0.000998 * q_dot_s_max_alt + 273.15;
elseif q_dot_s_max_alt > 7000 && q_dot_s_max_alt < 46000
    T = -23.4  - 0.00222 * q_dot_s_max_alt + 273.15;
else
    heights = [45.987, 47.948, 49.888, 51.825, 53.720, 55.646, 57.546] .* 1000;
    Temps = [148.78, 147.71, 146.72, 145.76, 144.75, 143.68, 142.68];
    T = interp1(heights, Temps, q_dot_s_max_alt);
end

Tupstream = T;
R = 8314/44.01;
gamma = 1.29;
M = q_dot_s_max_vel/sqrt(gamma * R * Tupstream)
Tdownstream = Tupstream * (((2*gamma*M^2-(gamma-1))*((gamma-1)*M^2+2))/((gamma+1)^2*(M^2)))
Tsurf = Tdownstream

% q_dot_rad = emissivity * stefanBoltzmannConst * Tsurf.^4
heatFluxIn = q_dot_s_max % - q_dot_rad
k = 1.63; % PICA - 1.63; SLA - 0.0576
Tvehicle = 150 + 273; % assume bondline temp of 150 degrees C
% thickness = -k * (Tvehicle - Tsurf) ./ heatFluxIn
% q_dot_s_max

% full energy method
rhoTPS = 0.265e3; %PICA - 0.265e3; SLA - 0.288e3 kg/m^3
CpTPS = 1625; %PICA - 1625; SLA - 1120 J/(kg*K)
maxT = 2200; %PICA - 2200; SLA - 2200 K
startingT = 200; % subject to change
thickness = max(cumtrapz(q_dot))/(rhoTPS*CpTPS*(maxT - startingT))

%% Local function
function [h_n_max, n_max, q_dot_s_max, J_s] = compute_performance(gamma, v_atm)


% EDL simulation constants
k = 1.9027e-4;
H = 11100; % m
mu_mars = 4.282837e13; % grav param of mars
R_mars = 3390e3; % radius of mars
r_orbiter = 370000 + R_mars;
rho0 = 0.02; % kg/m^3
g = 3.71; % m/s^2

% EDL design parameters
m = 100; % mass (kg)
r_b = 1; % radius of base (m)
A = pi*r_b; % reference area
r_n = 0.1; % radius of nosecone (m)
C_p_max = 2;
delta_c = 60 * (pi/180); % degrees - assume will change with more design
CD = C_p_max * (0.5 .* (1-sin(delta_c).^4)*(r_n./r_b).^2 + (sin(delta_c).^2).*(1-(r_n./r_b).^2.*cos(delta_c)));

beta = m/(CD*A); % ballistic coefficient
C = -(rho0.*H)./(2.*beta.*sin(gamma)); %

% Altitude of peak deceleration
h_n_max = H*log(-2*C);

% Magnitude of peak deceleration
n_max = ((v_atm.^2).*sin(gamma))./(53.3.*H);

% Peak heat rate
q_dot_s_max = (k ./ sqrt(r_n)) * ((beta .* sin(gamma) ./ (3 * H)).^0.5) .* ((v_atm.^3) ./ exp(0.5)); % (W/m^2)
% q_dot_s_max = q_dot_s_max./1e4; % W/cm^2

% Total integrated heat load
J_s = k .* v_atm.^2 .* ((beta*pi*H)./(rho0 .* r_n .* sin(gamma)))^0.5; % (J/m^2)
% J_s = J_s./1e4; % (J/cm^2)


end
