clear;clc;close all

t_span = 0:0.01:200;
h_initial = 1000;
v_initial = -30;
init_state = [h_initial; v_initial];
options = odeset('Events',@events);

states = []
thrusts = 3120;
t_delay = 5;
for i = 1:length(thrusts)
    [time, state] = ode45(@(t,x)eom(t,x,thrusts(i),t_delay(i)), t_span, init_state, options);

    states(i).state = state;
    states(i).time = time;
end

%%
figure(1)
for i = 1:length(thrusts)
    plot(states(i).time,states(i).state(:,1))
    hold on 
end
xlabel('Time (sec)')
ylabel('Altitude (m)')
%%
figure(2)
for i = 1:length(thrusts)
    plot(states(i).time,states(i).state(:,2))
    hold on 
end
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
%%
figure(5)
% for i = 1:length(chute_sizes)
%     plot(log(i).time,log(i).state(:,2), 'LineWidth',2)
%     hold on 
% end
% xlabel('\textbf{Time (sec)}','interpreter','latex','fontsize',12)
% ylabel('\textbf{Velocity (m/s)}','interpreter','latex','fontsize',12)
% title('\textbf{Nominal Chute Descent Velocity Profile}','interpreter','latex','fontsize',12)
% legend('5.0','7.5', '10.0','12.5', '15.0', '17.5', '20.0')


function [x_dot, log] = eom(t,x,T,t_delay)

    h = x(1);
    v = x(2);

    m = 606;
    g_mars = 3.73;

    if t < t_delay
        dhdt = v;
        dvdt = -g_mars;
        % dvdt = 0;
    elseif t > t_delay
        dhdt = v;
        dvdt = T/m - g_mars;
    end

    x_dot = [dhdt;dvdt];

    log.a = dvdt;

end

function [value,isterminal,direction] = events(t,x)
%rocket_events Called by ode45 to determine various edl events
value = x(1);
isterminal = 1;  % stop integration when event happens
direction = -1

end