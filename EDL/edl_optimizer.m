clear; clc; close all

% Initial guess
x0 = [5000; 1 * (pi/180)];

% lower and upper bounds
lb = [0, 20000];
ub = [0, 10 * (pi/180)];

options = optimoptions ('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                        'ConstraintTolerance', 1, ...
                        'StepTolerance', 1e-10);

% Solve the optimization problem
[ x_opt , fval ] = fmincon ( @edl_eom, x0 , [], [], [] , [] , lb,ub , @nonlinear_constraints4, options )

%% Local functions

function f = cost(x)
    
    v_atm = x(1); % re-entry velocity
    gamma = x(2); % re-entry angle

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
    m = 100; % mass (kg)
    r_b = 1; % radius of base (m)
    A = pi*r_b; % reference area
    v_atm = 5000; % reentry velocity (m/s)
    r_n = 0.1; % radius of nosecone (m)
    C_p_max = 2;
    delta_c = 60 * (pi/180); % degrees - assume will change with more design
    CD = C_p_max * (0.5 .* (1-sin(delta_c).^4)*(r_n/r_b)^2 + (sin(delta_c)^2)*(1-(r_n/r_b)^2*cos(delta_c)));
    
    beta = m/(CD*A); % ballistic coefficient
    C = -(rho0.*H)./(2.*beta.*sin(gamma)); %
    
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


    weight_n_max = 1;
    weight_q_dot_s_max = 0.0001;
    weight_J_s = 0.00001;
    f = weight_n_max * n_max + weight_q_dot_s_max * q_dot_s_max + weight_J_s * J_s;


end