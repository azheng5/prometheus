%% aerobraking.m
% Preliminary aerobraking simulations script for estimating time of flight
% and number of orbits of aerobraking phase. For higher fidelity
% aerobraking simulations, see aerobraking_high_fidel.py

clear;clc;close all

mu_m = 4.282837e4; % km^3/s^2
R_m = 3390; % radius of mars (km)

% Compute delta V's and velocity vectors in interplanetary transfer
interplanetary
v_sv2 = v2_in_km_s; % velocity of vehicle wrt sun at mars (km/s)
v_sp2 = vsp_Mars_km_s; % velocity of mars wrt sun (km/s)
v_inf2 = v_sv2 - v_sp2;
v_inf2_norm = norm(v_inf2)

h = 129; % altitude of aerobraking
r_p = h + R_m;
e_a1 = 0.9;

a_h1 = -mu_m/(v_inf2_norm^2); % a of arrival hyperbolic orbit (km)
a_a1 = r_p/(1-e_a1);
v_h1 = sqrt(mu_m*(2/r_p - 1/a_h1)); % velocity at periapsis of arrival hyperbolic orbit
v_a1 = sqrt(mu_m*(2/r_p - 1/a_a1)); % velocity at periapsis of initial aerobraking orbit
% v_a1 = 4; % from Lu
% v_a1 = sqrt(mu_m*(2/r_p - 1/r_p)); % going from MTO directly to LMO is 4.8 km/s
%overload
% v_h1 = 2.95
dv1 = v_h1 - v_a1 % delta-v from MTO to aerobraking orbit

% Calculate initial aerobraking orbit elements
incl_a1 = 70;
eps_a1 = (v_a1^2)/2 - mu_m/r_p;
a_a1 = -mu_m/(2*eps_a1);
e_a1 = (a_a1 - r_p)/a_a1;

% Calculate main aerobraking phase
num_aerobrake_orbits = 120; % ignore this value


% a_a1_meters = a_a1*1000;
% mu_m_meters = 4.282837e13;
% r_p_meters = r_p*1000
% a_aero = [];
% a_curr = a_a1_meters;
% % NOTE: these are done in meters not km
% for i = 1:num_aerobrake_orbits*(time_of_drag/dt)
% 
%     % meters
%     h = 129000;
%     T = -31  - 0.000998 * h;
%     p = 0.699 * exp(-0.00009*h);
%     rho = p / (0.1921 * (T + 273.15));
% 
%     v = sqrt(mu_m_meters*(2/r_p_meters - 1/a_curr));
%     a_dot = -2*(mu_m_meters^(-0.5))*(a_curr^1.5)*(rho*CD*A*(v^2)/(2*m));
%     a_next = a_curr + a_dot * dt;
% 
%     if mod(dt*i,time_of_drag) == 0
%         dt*i
%         a_aero = [a_aero; a_curr];
%     end
% 
%     if a_next <= 129000+R_m*1000
%         break
% 
%     end
%     a_curr= a_next;
% end
% a_final_meters = a_aero(end);
% a_final = a_final_meters/1000

dt = 60*5;
% T = 2*pi*sqrt((a_a1^3)/mu_m);
T = 42000;
% t_span = 1:dt:num_aerobrake_orbits*T;
t_span = [0 (0.155*365*24*60*60)];
init_state = [r_p;0;0;0;v_a1;0];
% options = odeset('RelTol', 0.1, 'AbsTol', 0.1);
[time, state] = ode89(@(t,x)two_body(t,x), t_span, init_state);

num_orbits = 0;
a_indices = [];
for i = 1:length(state)-1

    if state(i,2) < 0 && state(i+1,2) > 0
        num_orbits = num_orbits + 1;
        a_indices = [a_indices; i];
    end
    

end

orbit_span = 1:length(a_indices);
r_a_aero = zeros(length(a_indices),1);
for i = 1:length(a_indices)

    r_norm = norm(state(a_indices(i),1:3));
    v_norm = norm(state(a_indices(i),4:6));

    eps = (v_norm^2)/2 -mu_m/r_norm;
    a= -mu_m/(2*eps);
    e = (a - r_p)/a;
    r_a_aero(i) = a*(1+e);
    

    % if abs(r_a_aero(i) - r_p) < 1
    %     fprintf('required num of orbits: %d',i)
    % end

end

figure(1)
plot(time,state(:,1:3))
%%

hF = figure(2);
hA = axes(hF);

plot(state(:,1),state(:,2),'LineWidth',1,'Color','#9BBFF7')
xlabel('\textbf{x-position (km)}','interpreter','latex','fontsize',18)
ylabel('\textbf{y-position (km)}','interpreter','latex','fontsize',18)
title('\textbf{Aerobraking Sequence For 2 Months (100 Orbits) at 129 km altitude}','interpreter','latex','fontsize',18)
% set(hA, 'XTick', [], 'XTickLabel', []);
% set(hA, 'YTick', [], 'YTickLabel', []);
% set(get(hA, 'XAxis'), 'Visible', 'off');
% set(get(hA, 'YAxis'), 'Visible', 'off');
%%
figure(3)
plot(orbit_span,r_a_aero,'LineWidth',2)
hold on
plot(orbit_span,ones(length(orbit_span),1).*r_p,'LineWidth',2)
xlabel('\textbf{Number of Orbits}','interpreter','latex','fontsize',18)
ylabel('\textbf{Apoapsis}','interpreter','latex','fontsize',18)
title('\textbf{Aerobraking Sequence For 2 Months (100 Orbits) at 129 km altitude}','interpreter','latex','fontsize',18)
legend('Aerobraking Apoapsis','Reference Apoapsis')


function [x_dot] = two_body(t,x)
   t

   rx = x(1);
   ry = x(2);
   rz = x(3);
   vx = x(4);
   vy = x(5);
   vz = x(6);
   r = [rx;ry;rz];
   v = [vx;vy;vz];

   mu_m = 4.282837e4; % km^3/s^2
   R_m = 3390; % radius of mars (km)
   m = 1500; % rover mass
   CD = 1.28; % drag coefficient of a slab
   A = (5.492*1e-6)*2; % km^2
   beta = m/(CD*A);
   beta= beta/1e6; % 0.10669

   %%%%%


   h = norm(r) - R_m;
   %this section in meters%%%
   % h = 1290000
   h_meters = h*1000;
   % T = -23.4  - 0.00222 * h_meters;
   % p = 0.699 * exp(-0.00009*h_meters);
   % rho = p / (0.1921 * (T + 273.15));
   T = -31  - 0.000998 * h_meters;
   p = 0.699 * exp(-0.00009*h_meters);
   rho = p / (0.1921 * (T + 273.15));
   rho  = rho*1e9; % kg/km^3
   % rho = 0.02*1e9;

   % if h < 450
   %     heights = [45.987, 47.948, 49.888, 51.825, 53.720, 55.646, 57.546] .* 1000;
   %     densitys = [1.129e-4, 8.518e-5, 6.946e-5, 5.446e-5, 4.289e-5, 3.36e-5, 2.636e-05];
   %     rho = interp1(heights, densitys, h_meters,'linear','extrap')*1e9;
   % else
   %     rho = 0;
   % end


   F_D = -(0.5*CD*A*rho).*((norm(v))^2)*(v./norm(v));
   a_d = F_D./m;

   accel = (-mu_m .* r) ./ (norm(r)^3) + a_d;

   x_dot = [v;accel];


end