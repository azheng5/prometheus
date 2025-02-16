clear;clc;close all

% EDL simulation constants
k = 1.9027e-4;
H = 11100; % m
mu_mars = 4.282837e13; % grav param of mars
R_mars = 3390e3; % radius of mars
r_orbiter = 370000 + R_mars;
rho0 = 0.02; % kg/m^3
g = 3.71; % m/s^2

% EDL design parameters
gamma = [1:0.1:10] .* (pi/180); % flight path angle (rad)
gamma0 = 10 .* (pi/180);
m = 100; % mass (kg)
r_b = 1; % radius of base (m)
A = pi*r_b; % reference area
v_atm = 5000; % reentry velocity (m/s)
r_n = 0.1; % radius of nosecone (m)
C_p_max = 2;
delta_c = 60 * (pi/180); % degrees - assume will change with more design
CD = C_p_max * (0.5 .* (1-sin(delta_c).^4)*(r_n/r_b)^2 + (sin(delta_c)^2)*(1-(r_n/r_b)^2*cos(delta_c)));

beta = m/(CD*A); % ballistic coefficient
% beta = 90;
C = -(rho0.*H)./(2.*beta.*sin(gamma)); %

% Equations
% peak heat rate

% v_circ = sqrt(mu_mars/r_orbiter);
% v_deorbit = v_circ + delta_v;

% Altitude of peak deceleration
h_n_max = H*log(-2*C);

% Magnitude of peak deceleration
n_max = ((v_atm^2).*sin(gamma))./(53.3.*H);

% Peak heat rate
q_dot_s_max = (k / sqrt(r_n)) * ((beta .* sin(gamma) ./ (3 * H)).^0.5) .* ((v_atm^3) / exp(0.5)); % (W/m^2)
q_dot_s_max = q_dot_s_max/1e4; % W/cm^2

% Total integrated heat load
J_s = k .* v_atm^2 .* ((beta*pi*H)./(rho0 .* r_n .* sin(gamma))); % (J/m^2)
J_s = J_s/1e4; % (J/cm^2)

% Allen Eggers analytic solution
% v = v_atm * exp( ((-rho0/H)/(2*beta*sin(gamma))) * exp(-h/H) );

t_span = 0:0.1:5000;
init_state = [H;v_atm;gamma0];
config.beta = beta;
config.rho0 = rho0;
config.H = H;
config.R_mars = R_mars;
config.g = g;
options = odeset('Events',@edl_events);
[time, state] = ode45(@(t,x)edl_eom(t,x,config), t_span, init_state, options);

a = zeros(length(time),1);
for i = 1:length(time)
    [xdot,log] = edl_eom(time(i),state(i,:)',config);
    a(i) = log.a;
end



figure(1)
plot(gamma*(180/pi),q_dot_s_max)
xlabel('Flight path angle (deg)')
ylabel('Peak heating rate (W/cm^2)')

figure(2)
plot(gamma*(180/pi),J_s)
xlabel('Flight path angle (deg)')
ylabel('Total heat load (W/cm^2)')

figure(3)
plot(gamma*(180/pi),h_n_max)
xlabel('Flight path angle (deg)')
ylabel('Altitude of peak deceleration (m)')

figure(4)
plot(gamma*(180/pi),n_max)
xlabel('Flight path angle (deg)')
ylabel('Peak deceleration (m/s^2)')

figure(5)
plot(time,state(:,1))
xlabel('Time (sec)')
ylabel('Altitude (m)')

figure(6)
plot(time,state(:,2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')

figure(7)
plot(time,state(:,3) .* (180/pi))
xlabel('Time (sec)')
ylabel('Flight path angle (deg)')


% function rho = calc_rho(h,H)
% 
%     rho = rho0 * exp(-h/H);
% 
% end
