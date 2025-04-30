clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked')

%% Input Parameters

%Input Latitude
Phi = -76.5; %latitude (deg)

martianYear = 1:670;
Percent_sun = zeros(1,670);
rAU = zeros(1,670);
figure(1)

%Calling Function

for j = martianYear
    [Percent_sun(j),rAU(j)] = marsDaylight(j,Phi);
end

hold on
name = ['Lat = ',num2str(Phi),'^o'];
plot(martianYear,100*real(Percent_sun(:)),'LineWidth',1.2,'DisplayName',name)
%yline(45,'-',{'Minimum Daylit Percentage'})
xline(day2sol(2029,10,10),'-',{'Arrival Date: 10/15/29 '})
xline(day2sol(2030,08,09),'-',{'End of Lifetime: 08/19/30'})

Percent_sun_real = 100*real(Percent_sun);
csvwrite('SunlightPercentage.csv',[martianYear', Percent_sun_real'])

%Plotting
grid on
fontname('Times New Roman')
title('Sunlight Percentage over Martian Year')
xlabel('Date','FontWeight','bold')
ylabel('Daylit Percentage of Martian Sol','FontWeight','bold')

%Dates
sols = linspace(0, 670, 10); % Example sols for tick marks
start_date = datetime(2029, 3, 3); % Sol 0 corresponds to March 3, 2029
end_date = datetime(2031, 1, 19); % Sol 670 corresponds to January 19, 2031

% Convert sols to dates
days_per_sol = (end_date - start_date) / 670; % Find days per sol
dates = start_date + days_per_sol * sols; % Convert sols to corresponding dates

% Set custom xticks and labels
xticks(sols);
xticklabels(datestr(dates, 'mmm dd, yyyy')); % Format dates nicely

% Rotate labels for better readability
xtickangle(45);

%Solar Irradiance Plotting
Irradiance = 1366.1.*(1./rAU).^2;
figure(2)
plot(martianYear,Irradiance,'LineWidth',1.2)
grid on
title('Martian Solar Irradiance over a Martian Year')
xlabel('Martian Sol','FontWeight','bold')
ylabel('Solar Irradiance (W/m^2)','FontWeight','bold')
xline(day2sol(2029,10,10),'-',{'Arrival Date'})
xline(day2sol(2030,08,09),'-',{'Mission End Date'})


xticks(sols);
xticklabels(datestr(dates, 'mmm dd, yyyy')); % Format dates nicely

% Rotate labels for better readability
xtickangle(45);


%% Function Definition

function [L_sun,rAU] = marsDaylight(sol_num,Phi)
%Martian Parameters
mu_s = 1.327 * 10^11; %km^3/s^2
a = 227936640; %km
e = 0.0935;
f_w = deg2rad(170.25); %solstice anomaly
theta_a = -deg2rad(24.936); %Martian axial tilt

%Time
%Martian Day 0 is defined as Northern Summer Solstice, which is also
%Southern Winter Solstice

%Enter Martian Day (1 Martian Year is 670 Martian Sols)
martianDayLength = 88775.245; %sec
t = sol_num * martianDayLength; %sec

%Tau = time at perihelion
tau = 0; %sec

%Mean anomaly
n = sqrt(mu_s/a^3); %s^-1
M = n*(t-tau); %radians

%True anomaly
f = M + (2*e - 0.25*e^3)*sin(M) + 1.25*e^2*sin(2*M)+(13/12)*(e^3*sin(3*M)); %rad

% Angular Momentum
h = sqrt(mu_s * a*(1-e^2)); %km^2/s

%Radius of Orbit in AU
r = a*(1-e^2)/(1 + e*cos(f)); %km
rAU = r/149597870.7;

%Solar declination angle
delta = rad2deg(theta_a * cos(f-f_w)); %deg

%Day Length Angle
w0 = 2*acosd(-tand(Phi)*tand(delta)); %deg

%Daylight Length in % of Day Length
L_sun = w0/360;

end
