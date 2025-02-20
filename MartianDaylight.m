clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked')

%% Input Parameters

%Input Latitude
Phi = -65:-5:-80; %latitude (deg)

martianYear = 1:687;
Percent_sun = zeros(1,687);
figure(1)

%Calling Function
for i = 1:length(Phi)
    for j = martianYear
        Percent_sun(i,j) = marsDaylight(j,Phi(i));
    end

    hold on
    name = ['Lat = ',num2str(Phi(i)),'^o'];
    plot(martianYear,100*real(Percent_sun(i,:)),'LineWidth',1.2,'DisplayName',name)
end


%Plotting
grid on
title('Sunlight Percentage over Martian Year')
subtitle('Day 0 is Martian Northern Summer Solstice')
xlabel('Martian Sol','FontWeight','bold')
ylabel('Daylit Percentage of Martian Sol','FontWeight','bold')
legend show


%% Function Definition

function L_sun = marsDaylight(sol_num,Phi)
%Martian Parameters
mu_s = 1.327 * 10^11; %km^3/s^2
a = 227936640; %km
e = 0.0935;
f_w = deg2rad(0); %solstice anomaly
theta_a = -deg2rad(24.936); %Martian axial tilt

%Time
%Martian Day 0 is defined as Northern Summer Solstice, which is also
%Southern Winter Solstice

%Enter Martian Day (1 Martian Year is 687 Martian Sols)
martianDayLength = 88775.245; %sec
t = sol_num * martianDayLength; %sec

%Tau = time at perihelion
tau = -24857068.35; %sec

%Mean anomaly
n = sqrt(mu_s/a^3); %s^-1
M = n*(t-tau); %radians

%True anomaly
f = M + (2*e - 0.25*e^3)*sin(M) + 1.25*e^2*sin(2*M)+(13/12)*(e^3*sin(3*M)); %rad

% Angular Momentum
h = sqrt(mu_s * a*(1-e^2)); %km^2/s

%Radius of Orbit
r = a*(1-e^2)/(1 + e*cos(f)); %km

%Solar declination angle
delta = rad2deg(theta_a * cos(f-f_w)); %deg

%Day Length Angle
w0 = 2*acosd(-tand(Phi)*tand(delta)); %deg

%Daylight Length in % of Day Length
L_sun = ((w0/360)*martianDayLength)/martianDayLength;

end
