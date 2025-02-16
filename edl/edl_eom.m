function [x_dot, log] = edl_eom(t,x,config)

    beta = config.beta;
    H = config.H;
    R_mars = config.R_mars;
    g = config.g;
    rho0 = config.rho0;

    h = x(1);
    v = x(2);
    gamma = x(3);

    rho = rho0 * exp(-h/H);

    dhdt = -v*sin(gamma);
    dvdt = -(rho * v^2)/(2*beta) + g*sin(gamma);
%     dgammadt = -(v * cos(gamma))/(R_mars + h) + g*cos(gamma)/v;
    dgammadt = 0;

    x_dot = [dhdt;dvdt;dgammadt];

    log.a = dvdt;


end