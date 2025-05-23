% Script that computes deorbit delta-V necessary to achieve a desired
% re-entry velocity into Mars

clear;clc;close all

mu_m = 4.282837e13; % grav param of mars (m^3/s^2)
R_m = 3390e3; % radius of mars
altitude = 129000;
r_circ = altitude + R_m;
r_op = 370000 + R_m;
h = 57000;

v_circ = sqrt(mu_m/r_circ);

% v_entry = [1439:1:3600];
v_entry = [1000:1:3540];
r_entry = R_m + h; % re-entry radius
eps_entry = 0.5.*v_entry.^2 - mu_m ./ r_entry;
a_entry = -mu_m./(2.*eps_entry);
e_entry = (r_circ - a_entry)./a_entry;
rp_entry = a_entry .* (1 - e_entry);
ra_entry = a_entry .* (1 + e_entry);
f_entry = acosd( (a_entry.*(1-e_entry.^2) - r_entry) ./ (r_entry.*e_entry)  );
% giving out of [-1 1] for some reason -> bc some re entry velocities are
% not feasible
% f_deorbit = acos( (a_entry*(1-e_entry^2) - r_circ) / (r_circ*e_entry)  ); % should be pi
v_deorbit = sqrt(mu_m .* (2./r_circ - 1./a_entry));
gamma = atan2d(e_entry .* sind(f_entry), 1 + e_entry .* cosd(f_entry));

delta_v = abs(v_deorbit - v_circ);

% Compute orbit raise delta v
deorbit_delta_v = 275;
a_t = (r_circ + r_op)/2;
delta_v_1 =  sqrt(mu_m*(2/r_circ - 1/a_t)) - sqrt(mu_m/r_circ);
delta_v_2 = sqrt(mu_m/r_op) - sqrt(mu_m*(2/r_op - 1/a_t));
orbit_raise_delta_v = delta_v_1 + delta_v_2;


%%
figure(1)
plot(v_entry, f_entry)
xlabel('Entry velocity (m/s)')
ylabel('True anomaly at re-entry point (deg)')
title('Feasible re-entry velocities from a 370 km parking orbit')
grid on

%%
figure(2)
plot(v_entry, delta_v,'Color','k','LineWidth',2)
xlabel('\textbf{Entry Velocity (m/s)}','interpreter','latex','fontsize',12)
ylabel('\textbf{Required} $\mathbf{\Delta V}$','interpreter','latex','fontsize',12)
title('\textbf{Deorbit Burn for Different Re-entry Velocities}','interpreter','latex','fontsize',12)
% grid on
% ax = gca;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
% ax.MinorGridLineStyle = ':';
% ax.XMinorGrid = 'on';
% ax.YMinorGrid = 'on';

%%
figure(3)
plot(v_entry, gamma,'Color','k','LineWidth',2)
xlabel('\textbf{Entry Velocity (m/s)}','interpreter','latex','fontsize',12)
ylabel('\textbf{Resultant Flight Path Angle (deg)}','interpreter','latex','fontsize',12)
title('\textbf{Feasible Re-entry Parameters From Orbit Geometry}','interpreter','latex','fontsize',12)
% grid on
% ax = gca;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
% ax.MinorGridLineStyle = ':';
% ax.XMinorGrid = 'on';
% ax.YMinorGrid = 'on';
%%
figure(4)
plot(f_entry, gamma)
ylabel('Flight path angle (deg)')
xlabel('True anomaly at re-entry point (deg)')
title('Feasible re-entry velocities from a 370 km parking orbit')
grid on

% %%
% figure(4)
% plot(delta_v, a_entry)
% xlabel('Delta-v burn (m/s)')
% ylabel('Entry trajectory semi-major axis (m)')
% title('Effect of delta v on semi-major axis')
% grid on
% 
% %%
% figure(5)
% plot(delta_v, ra_entry)
% xlabel('Delta-v burn (m/s)')
% ylabel('Entry trajectory apoapsis (m)')
% title('Effect of delta v on apoapsis')
% grid on
% 
% %%
% figure(6)
% plot(delta_v, rp_entry)
% xlabel('Delta-v burn (m/s)')
% ylabel('Entry trajectory periapsis (m)')
% title('Effect of delta v on periapsis')
% grid on