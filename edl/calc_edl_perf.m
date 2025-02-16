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
m = 200; % mass (kg)
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
J_s = k .* v_atm.^2 .* ((beta*pi*H)./(rho0 .* r_n .* sin(gamma))); % (J/m^2)
% J_s = J_s./1e4; % (J/cm^2)


end