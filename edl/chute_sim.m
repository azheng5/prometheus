clear;clc;close all

t_span = 0:0.01:2000;
h_initial = 15700;
v_initial = 500;
init_state = [h_initial; v_initial];
options = odeset('Events',@chute_events);
% chute_sizes = [5 7.5 10 12.5 15 17.5 20];
chute_sizes = 20;
log(length(chute_sizes)) = struct();
states = zeros(length(t_span),2*length(chute_sizes));
for i = 1:length(chute_sizes)
    [time, state] = ode45(@(t,x)chute_eom(t,x,chute_sizes(i)), t_span, init_state, options);
    chute_sizes(i)
    h_final = state(end,1)
    v_final = state(end,2)

    log(i).time = time;
    log(i).state = state;
end

figure(4)
for i = 1:length(chute_sizes)
    plot(log(i).time,log(i).state(:,1))
    hold on 
end
xlabel('Time (sec)')
ylabel('Altitude (m)')
% legend('5.0','7.5', '10.0','12.5', '15.0', '17.5', '20.0')

figure(5)
for i = 1:length(chute_sizes)
    plot(log(i).time,log(i).state(:,2), 'LineWidth',2)
    hold on 
end
xlabel('\textbf{Time (sec)}','interpreter','latex','fontsize',12)
ylabel('\textbf{Velocity (m/s)}','interpreter','latex','fontsize',12)
title('\textbf{Nominal Chute Descent Velocity Profile}','interpreter','latex','fontsize',12)
% legend('5.0','7.5', '10.0','12.5', '15.0', '17.5', '20.0')

%% For concatenating chute sim and ballistic descent sim
cat_time = [time; log.time+time(end)];
cat_state = [state(:,1:2);log.state];

% figure(1)
% plot(cat_time,cat_state(:,1),'LineWidth',2,'Color','b')
% xlabel('\textbf{Time (sec)}','interpreter','latex','fontsize',12)
% yyaxis left
% ylabel('\textbf{Altitude (m)}','interpreter','latex','fontsize',12)
% ax = gca;
% ax.YAxis(1).Color = 'b';
% 
% yyaxis right
% plot(cat_time,cat_state(:,2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
% ylabel('\textbf{Velocity (m/s)}','interpreter','latex','fontsize',12)
% ax = gca;
% ax.YAxis(2).Color = [0.8500 0.3250 0.0980];
% title('\textbf{Endurance EDL Profile}','interpreter','latex','fontsize',12)
%%

figure(1)
subplot(2,1,1)
plot(cat_time,cat_state(:,1),'LineWidth',2,'Color','k')
% xlabel('\textbf{Time (sec)}','interpreter','latex','fontsize',12)
% ylabel('\textbf{Altitude (m)}','interpreter','latex','fontsize',12)
ax = gca;

subplot(2,1,2)
plot(cat_time,cat_state(:,2),'LineWidth',2,'Color','k')
% xlabel('\textbf{Time (sec)}','interpreter','latex','fontsize',12)
% ylabel('\textbf{Velocity (m/s)}','interpreter','latex','fontsize',12)
% title('\textbf{Endurance EDL Velocity Profile}','interpreter','latex','fontsize',12)

%%
function [x_dot] = chute_eom(t,x,d)

    h = x(1);
    v = x(2);

    H = 11100; % m
    max_q = 750; % Pa
    rho0 = 0.02; % kg/m^3
    rho = rho0 * exp(-h/H);
    m = 900;
    CD = 1.2;
    g_mars = 3.73;
    A = pi*(d/2)^2;
    F_D = 0.5 * rho * (v^2) * CD * A;

    if h > 50
        dhdt = -v;
        dvdt = (-1/m)*F_D + g_mars;
    else
        dhdt = -v;
        dvdt = g_mars;
    end


    x_dot = [dhdt;dvdt];

end

function [value,isterminal,direction] = chute_events(t,x)
%rocket_events Called by ode45 to determine various edl events
    % i.e. capsule hitting ground
% value = x(1)-50; % define event as h = 0
value = x(1)-1000;
isterminal = 1;  % stop integration when event happens
direction = -1;   % only stop if h is decreasing

end