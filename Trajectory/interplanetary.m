%% interplanetary.m
% Script for calculating total delta V required to conduct interplanetary
% transfer from Earth to Mars in 2028-2029 transfer period. Relies on a
% Gauss' problem solver. Note: Running the simulation has a dependency on
% mice and the associated ephemerides specified in the meta-kernel file
% ephemerides.tm.

% Reset workspace
% clear; clc; close all
cspice_kclear

% Read metakernel file
mkfile = 'ephemerides.tm';
cspice_furnsh(mkfile); % load kernel file into matlab

% Create time
t0 = cspice_str2et('2028 DEC 22 00:00:00.00');% departure from Earth SOI
tf = cspice_str2et('2029 AUG 10 00:00:00.00'); % arrival to Mars SOI

% Get earth (e) with respect to solar system barycenter (m) in ECLIPJ2000 frame
[x_ES1, ~] = cspice_spkezr('EARTH', t0, 'ECLIPJ2000', 'NONE', '0');
p_ES1 = x_ES1(1:3);
v_ES1 = x_ES1(4:6);

% % Get mars (e) with respect to solar system barycenter (m) in ECLIPJ2000 frame
% [x_MS1, ~] = cspice_spkezr('MARS', t0, 'ECLIPJ2000', 'NONE', '0');
% p_MS1 = x_MS1(1:3);
% v_MS1 = x_MS1(4:6);

% Get mars (e) with respect to solar system barycenter (m) in ECLIPJ2000 frame
[x_MS2, ~] = cspice_spkezr('4', tf, 'ECLIPJ2000', 'NONE', '0');
p_MS2 = x_MS2(1:3);
v_MS2 = x_MS2(4:6);

% r_Earth = 1; % AU
% r_Mars = 1.524; % AU

% dir_Earth = [1; 0; 0];
% dir_Mars = [-cos(20*(pi/180)); -sin(20*(pi/180)); 0];

% r1 = r_Earth .* dir_Earth;
% r2 = r_Mars .* dir_Mars;


% Assumptions: coplanar Earth/Mars, 7.6 months transfer time, 200-deg
% transfer angle (as compared to a 180 degree Hohmann transfer angle)
%% Setup Gauss' problem
r1 = p_ES1 ./ 1.496e8; % convert from km to AU
r2 = p_MS2 ./ 1.496e8;
isShortWay = true;
alpha = 0.01;
t_in_seconds = tf - t0; %60*60*24*30*8.5;
t_in_TU = (1.9909e-7)*(t_in_seconds);

[v1, v2] = Gauss_p_iteration(r1,r2,t_in_TU,isShortWay,alpha)

convert_canonical_to_km_s = (149.598e6)/(5.0229e6);

v1_in_km_s = v1*(convert_canonical_to_km_s);
v2_in_km_s = v2*(convert_canonical_to_km_s)

% vsp_dir = [sin(20*(pi/180)); -cos(20*(pi/180)); 0];

mu_Sun = 1.327e11;
% r_Mars_km = r_Mars*149.598e6;
% vsp_mag_Mars = sqrt(mu_Sun/r_Mars_km);

% vsp_Mars_km_s = vsp_mag_Mars.*vsp_dir
vsp_Mars_km_s = v_MS2;
vsp_mag_Mars_km_s = norm(v_MS2);

v_inf_Mars_km_s = norm(vsp_Mars_km_s - v2_in_km_s)

R_e = 6378; % km
mu_e = 3.986e5; %km^3/s^2
r_p = 200 + R_e; % km
v_inf1 = norm(v_ES1 - v1_in_km_s);
a1 = -mu_e/(v_inf1^2);
v_p1 = sqrt(mu_e*(2/r_p - 1/a1));
v_c1 = sqrt(mu_e/r_p);
delta_v_1 = v_p1 - v_c1;


% Assume circular Martian parking orbit at the desired inclination

% mu_Mars = 4.282837e4; %km^3/s^2
% r_parking_orbit = 3396.2 + 370; % km
% v_parking_orbit = sqrt(mu_Mars/(3396.2 + 370));

% DV_Parking_Injection = v_inf_Mars_km_s - v_parking_orbit % km/s


% r_Mars_surface = 3396.2; % km
% v_Mars_surface_orbit = sqrt(mu_Mars/r_Mars_surface)

% a = (r_Mars_surface+r_parking_orbit)/2;
% v_transfer = sqrt(mu_Mars*(2/r_parking_orbit - 1/a))
% DV_Disposal = v_parking_orbit-v_transfer % km/s

% Now, just need to quantify our delta V to enter and exit the de-orbit
% burn for EDL stuff.


function[v1, v2] = Gauss_p_iteration(r1, r2, t, isShortWay, alpha)

%% Inputs: r1 and r2 vectors in canonical variables, t in canonical variables, alpha (regulator for Newton-Raphson)
%% Outputs: v1 and v2 in DU/TU canonical variables

r1 = r1;
r2 = r2;
t = t;

mu = 1; % due to usage of canonical variables
cosdf = dot(r1, r2)/(norm(r1)*norm(r2));
delta_f = acosd(cosdf);
delta_f = delta_f*pi/180;

if ~isShortWay
    delta_f = 2*pi - delta_f;
end

k = norm(r1)*norm(r2)*(1-cos(delta_f));
l = norm(r1) + norm(r2);
m = norm(r1)*norm(r2)*(1 + cos(delta_f));

F_func = @(p, delta_f, r2) 1 - (r2/p)*(1-cos(delta_f)); %Pass in norm(r2)
Fdot_func = @(p, mu, delta_f, r1, r2) sqrt(mu/p)*tan(delta_f/2)*((1-cos(delta_f))/p - 1/r1 - 1/r2); %Pass in norm(r)'s
G_func = @(p, r1, r2, delta_f, mu) (r1*r2*sin(delta_f))/sqrt(mu*p);

p_guess = ((k/(l - sqrt(2*m))) + (k/(l + sqrt(2*m))))/2;

p_k = p_guess;
a_func = @(p, m, l, k) (m*k*p)/((2*m - l^2)*p^2 + 2*k*l*p - k^2);
a_k = a_func(p_k, m, l, k);

if a_k > 0

    deltaE_func = @(p, delta_f, r1, r2, mu, a, F_func, Fdot_func) atan2((-r1*r2*Fdot_func(p, mu, delta_f, r1, r2))/sqrt(mu*a), 1 - (r1/a)*(1-F_func(p, delta_f, r2)));
    deltaE = wrapTo2Pi(deltaE_func(p_k, delta_f, norm(r1), norm(r2), mu, a_k, F_func, Fdot_func));
    t_func = @(p, G_func, a, mu, deltaE, r1, r2, delta_f) G_func(p, r1, r2, delta_f, mu) + sqrt((a^3)/mu)*(deltaE - sin(deltaE));
    dt_dp_func = @(p, G_func, r1, r2, delta_f, mu, a, t, k, m, l, deltaE) -G_func(p, r1, r2, delta_f, mu)/(2*p) - 1.5*a*(t-G_func(p, r1, r2, delta_f, mu))*((k^2 + (2*m - l^2)*(p^2))/(m*k*(p^2))) + sqrt((a^3)/mu)*((2*k*sin(deltaE))/(p*(k-l*p)));
    t_k = t_func(p_k, G_func, a_k, mu, deltaE, norm(r1), norm(r2), delta_f);
    dt_dp_k = dt_dp_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, t_k, k, m, l, deltaE);

else % a_k < 0

    deltaF_func = @(p, r1, a, F_func, r2, delta_f) acosh(1 - (r1/a)*(1-F_func(p, delta_f, r2)));
    deltaF = deltaF_func(p_k, norm(r1), a_k, F_func, norm(r2), delta_f);
    t_func = @(p, G_func, r1, r2, delta_f, mu, a, deltaF) G_func(p, norm(r1), norm(r2), delta_f, mu) + sqrt(((-a)^3)/mu)*(sinh(deltaF) - deltaF);
    dt_dp_func = @(p, G_func, r1, r2, delta_f, mu, a, t, k, m, l, deltaF) -G_func(p, r1, r2, delta_f, mu)/(2*p) - 1.5*a*(t-G_func(p, r1, r2, delta_f, mu))*((k^2 + (2*m - l^2)*(p^2))/(m*k*(p^2))) - sqrt(((-a)^3)/mu)*((2*k*sinh(deltaF))/(p*(k-l*p)));
    t_k = t_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, deltaF);
    dt_dp_k = dt_dp_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, t_k, k, m, l, deltaF);
end

tol = 10e-6;
alpha = alpha;
while abs(t - t_k) >= tol

    p_next = p_k + alpha*((t - t_k)/dt_dp_k);

    p_k = p_next;
    a_k = a_func(p_k, m, l, k);

    if a_k > 0
    
        deltaE_func = @(p, delta_f, r1, r2, mu, a, F_func, Fdot_func) atan2((-r1*r2*Fdot_func(p, mu, delta_f, r1, r2))/sqrt(mu*a), 1 - (r1/a)*(1-F_func(p, delta_f, r2)));
        deltaE = wrapTo2Pi(deltaE_func(p_k, delta_f, norm(r1), norm(r2), mu, a_k, F_func, Fdot_func));
        t_func = @(p, G_func, a, mu, deltaE, r1, r2, delta_f) G_func(p, r1, r2, delta_f, mu) + sqrt((a^3)/mu)*(deltaE - sin(deltaE));
        dt_dp_func = @(p, G_func, r1, r2, delta_f, mu, a, t, k, m, l, deltaE) -G_func(p, r1, r2, delta_f, mu)/(2*p) - 1.5*a*(t-G_func(p, r1, r2, delta_f, mu))*((k^2 + (2*m - l^2)*(p^2))/(m*k*(p^2))) + sqrt((a^3)/mu)*((2*k*sin(deltaE))/(p*(k-l*p)));
        t_k = t_func(p_k, G_func, a_k, mu, deltaE, norm(r1), norm(r2), delta_f);
        dt_dp_k = dt_dp_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, t_k, k, m, l, deltaE);

    else % a_k < 0

        deltaF_func = @(p, r1, a, F_func, r2, delta_f) acosh(1 - (r1/a)*(1-F_func(p, delta_f, r2)));
        deltaF = deltaF_func(p_k, norm(r1), a_k, F_func, norm(r2), delta_f);
        t_func = @(p, G_func, r1, r2, delta_f, mu, a, deltaF) G_func(p, norm(r1), norm(r2), delta_f, mu) + sqrt(((-a)^3)/mu)*(sinh(deltaF) - deltaF);
        dt_dp_func = @(p, G_func, r1, r2, delta_f, mu, a, t, k, m, l, deltaF) -G_func(p, r1, r2, delta_f, mu)/(2*p) - 1.5*a*(t-G_func(p, r1, r2, delta_f, mu))*((k^2 + (2*m - l^2)*(p^2))/(m*k*(p^2))) - sqrt(((-a)^3)/mu)*((2*k*sinh(deltaF))/(p*(k-l*p)));
        t_k = t_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, deltaF);
        dt_dp_k = dt_dp_func(p_k, G_func, norm(r1), norm(r2), delta_f, mu, a_k, t_k, k, m, l, deltaF);
    end

end

t = t_k;
p = p_k;

F = F_func(p, delta_f, norm(r2));
Fdot = Fdot_func(p, mu, delta_f, norm(r1), norm(r2));
G = G_func(p, norm(r1), norm(r2), delta_f, mu);
Gdot = (1 + Fdot*G)/F;

v1 = (r2 - F*r1)/G;
v2 = (Gdot*r2 - r1)/G;

end