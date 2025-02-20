clear;clc;close all
v_atm = [1600:1:3600]; % reentry velocity (m/s)
[h_n_max, n_max, q_dot_s_max, J_s] = compute_performance(5*(pi/180), v_atm);

emissivity = 0.9; % PICA - 0.9; SLA - 0.7
stefanBoltzmannConst = 5.67e-8;

Tsurf = 1600 + 273;
% q_dot_rad = emissivity * stefanBoltzmannConst * Tsurf.^4
heatFluxIn = q_dot_s_max*10^4 % - q_dot_rad
k = 1.63; % PICA - 1.63; SLA - 0.0576
Tvehicle = 150 + 273; % assume bondline temp of 150 degrees C
thickness = -k * (Tvehicle - Tsurf) ./ heatFluxIn

figure(1)
plot(v_atm,q_dot_s_max)
xlabel('Re-entry veloicty (m/s)')
ylabel('Peak heating rate (W/cm^2)')

figure(2)
plot(v_atm,J_s)
xlabel('Re-entry veloicty (m/s)')
ylabel('Total heat load (W/cm^2)')

figure(3)
plot(v_atm,h_n_max)
xlabel('Re-entry veloicty (m/s)')
ylabel('Altitude of peak deceleration (m)')

figure(4)
plot(v_atm,n_max)
xlabel('Re-entry veloicty (m/s)')
ylabel('Peak deceleration (m/s^2)')

figure(5)
plot(v_atm,thickness)
xlabel('Re-entry velocity (m/s)')
ylabel('Heat shield thickness (m)')

%%

gamma = [1:0.1:80] .* (pi/180); % flight path angle (rad)
v_atm = 5000;

[h_n_max, n_max, q_dot_s_max, J_s] = compute_performance(gamma, v_atm);


figure(6)
plot(gamma * (180/pi),q_dot_s_max)
xlabel('Flight path angle (deg)')
ylabel('Peak heating rate (W/cm^2)')

figure(7)
plot(gamma * (180/pi),J_s)
xlabel('Flight path angle (deg)')
ylabel('Total heat load (W/cm^2)')

figure(8)
plot(gamma * (180/pi),h_n_max)
xlabel('Flight path angle (deg)')
ylabel('Altitude of peak deceleration (m)')

figure(9)
plot(gamma * (180/pi),n_max)
xlabel('Flight path angle (deg)')
ylabel('Peak deceleration (m/s^2)')


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
q_dot_s_max = q_dot_s_max./1e4; % W/cm^2

% Total integrated heat load
J_s = k .* v_atm.^2 .* ((beta*pi*H)./(rho0 .* r_n .* sin(gamma))); % (J/m^2)
J_s = J_s./1e4; % (J/cm^2)


end