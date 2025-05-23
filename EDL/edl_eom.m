function [x_dot, log] = edl_eom(t,x,config)

    beta = config.beta;
    H = config.H;
    R_mars = config.R_mars;
    g = config.g;
    rho0 = config.rho0;
    k = config.k;
    r_n = config.r_n;
    m = config.m;

    h = x(1);
    v = x(2);
    gamma = x(3);
    q_s = x(4);


    if h < 7000
        T = -31  - 0.000998 * h;
        p = 0.699 * exp(-0.00009*h);
        rho = p / (0.1921 * (T + 273.15));
    elseif h > 7000%% &&  h < 46000
        T = -23.4  - 0.00222 * h;
        p = 0.699 * exp(-0.00009*h);
        rho = p / (0.1921 * (T + 273.15));

    else

        heights = [45.987, 47.948, 49.888, 51.825, 53.720, 55.646, 57.546] .* 1000;
        %     p = [3.2, 2.490E-01, 1.510E-01, 1.050E-01, 9.170E-02, 7.140E-02];
        %     t = [147.71, 146.25, 145.76, 143.68, 142.68];
        densitys = [1.129e-4, 8.518e-5, 6.946e-5, 5.446e-5, 4.289e-5, 3.36e-5, 2.636e-05];
        rho = interp1(heights, densitys, h);

    end

%     max_q = 750; % Pa

%     rho = rho0 * exp(-h/H);


    dhdt = -v*sin(gamma);
    dvdt = -(rho * v^2)/(2*beta) + g*sin(gamma);
    dgammadt = (1/v) * (-((v^2) * cos(gamma))/(R_mars + h) + g*cos(gamma) );
%     dgammadt = 0;
    dqdt = k * ((rho/r_n)^0.5) * (v^3);

%     C_D_chute = 0.65;
%     d_chute = 30;
%     rho = 0.006372042676251;
%     F_D = 0.5 * rho * (v^2) * C_D_chute * pi* (d_chute/2)^2;
%     sqrt(max_q/(0.5*rho))
%     if v < sqrt(max_q/(0.5*rho)) && h < 7000
%         dhdt = -v;
%         dvdt = -(1/m)*F_D;
%     end

    x_dot = [dhdt;dvdt;dgammadt;dqdt];

    dynamic_p = 0.5*rho*(v^2);

    log.a = dvdt;
    log.q_dot = dqdt;
    log.dynamic_p = dynamic_p;
    log.rho = rho;


end