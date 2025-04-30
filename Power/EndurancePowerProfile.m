clear; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

mDayLength = 24.661*60; %min

%Duration of Each mode
rovSci = mDayLength * 0.16062;
rovSam = mDayLength * 0.12165;
rovMov = mDayLength * 0.08110;
rovSaf = mDayLength * 0.05330;
%rovCal = mDayLength * 0.03590;
rovEcl = mDayLength * 0.58333;

%Beginning times of each mode
timeSci = 0; %min
timeMov = timeSci + rovSci; %min
timeSam = timeMov + rovMov; %min
%timeCal = timeSaf + rovSaf;
timeSaf = timeSam + rovSam;
timeEcl = timeSaf + rovSaf; %min
timeEnd = timeEcl + rovEcl; %min

%Power consumption of each mode
pSci = 144.24548; %W
pSam = 143.18948; %W
pMov = 205.9972493; %W
pSaf = 129.3728; %W
%pCal = 166.8984; %W
pEcl = 165.6715328; %W
pEnd = pSci; %W

time = [timeSci,timeMov,timeSam,timeSaf,timeEcl,timeEnd];
power = [pSci,pMov,pSam,pSaf,pEcl,pEcl];

figure(1)
stairs(time/60,power,'LineWidth',1.2)
grid off
fontname('Times New Roman')
xlabel('Time (hrs)','FontWeight','bold','fontsize',12,'interpreter','latex')
ylabel('Power Draw (W)','FontWeight','bold','fontsize',12,'interpreter','latex')
xlim([0,timeEnd/60])
ylim([120,220])
title('Endurance Power Draw per Martian Sol','fontsize',12,'interpreter','latex')
% xline(timeSci/60,'Label','Science Mode','LabelVerticalAlignment','top')
% xline(timeSam/60,'Label','Sample Mode','LabelVerticalAlignment','top')
% xline(timeMov/60,'Label','Movement Mode','LabelVerticalAlignment','top')
% xline(timeSaf/60,'Label','Safe/Charge Mode','LabelVerticalAlignment','top')
% xline(timeCal/60,'Label','Calibration Mode','LabelVerticalAlignment','top')
% xline(timeEcl/60,'Label','Eclipse Mode','LabelVerticalAlignment','top')
