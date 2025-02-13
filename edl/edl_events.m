function [value,isterminal,direction] = edl_events(t,x)
%rocket_events Called by ode45 to determine various edl events
    % i.e. capsule hitting ground
value = x(1); % define event as h = 0
isterminal = 1;  % stop integration when event happens
direction = -1;   % only stop if h is decreasing

end