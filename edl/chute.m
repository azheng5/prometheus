% Computes parachute sizing

clear;clc; close all

h_descent = 10000; % rough height of chute descent phase
T = -31  - 0.000998 * h_descent;
p = 0.699 * exp(-0.00009*h_descent);
rho = p / (0.1921 * (T + 273.15));
m = 900; % rover mass

% Requirments
max_M = 2.3; % speed of sound 240 m/s
max_q = 750; % Pa
% q_chute = 0.5*rho*v_initial^2 % (Pa) dynamic pressure at chute deployment

v_initial = sqrt(max_q/(0.5*rho)) % 485.184(m/s) velocity at moment of chute deployment
h_initial =  6414; % simulated altitude at v_initial of 485.184 m/s for 1500 kg capsule. 

% sort of arbitrary, but we want chutes on as long as possible without
% imposing a risk that the chute will fall on the rover. MSL ended chute
% descent at 1500 m but it used powered descent, so the altitude should be
% quite a bit lower than that
% if ballistic descent continued and chute didn't deploy at time it was 
% suposed to, rover would hit ground 30 sec later.
h_final = 50;

% Desired velocity to slow chute down (at h_final)
% Assuming drop altitude of 100 m, 
g_mars = 3.73;

max_v_impact = 25; % slightly higher than Spirit
t_freefall = 5; % min time to ensure chute doesn't hit rover
delta_h = h_initial - h_final; % (m) distance of parachute descent phase
% v_desired = max_v_impact - g_mars*t_freefall
v_desired = 10.167;
t_freefall = (max_v_impact - v_desired)/g_mars;

% height 
% t_freefall = (max_v_impact - v_desired)/g_mars;
% CD = 0.85; % CD of parachute

% v_avg = (v_initial - v_desired)/2;
% delta_t = delta_h/v_avg;
% F_D = (m*(v_initial - v_desired))/delta_t;

% d_chute = sqrt((8*F_D)/(pi*rho*(v_avg^2)*CD)) % (m) diameter of chute