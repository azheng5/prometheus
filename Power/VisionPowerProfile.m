clear; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

orbPeriod = 7017.28; %sec

%Duration of Each mode
orbNom = orbPeriod * 0.609 + 0.03;
%orbCal = orbPeriod * 0.03;
orbEcl = orbPeriod * 0.358;

%Beginning times of each mode
timeNom = 0; %sec
%timeCal = timeNom + orbNom; %sec
timeEcl = timeNom + orbNom; %sec
timeEnd = timeEcl + orbEcl;

%Power consumption of each mode
pNom = 1357.402533; %W
%pCal = 1106.078638; %W
pEcl = 561.3623765; %W
pEnd = pNom;

time = [timeNom,timeEcl,timeEnd];
power = [pNom, pEcl, pEcl];

figure(1)
stairs(time/60,power,'LineWidth',1.2)
grid off
fontname('Times New Roman')
xlabel('Time (min)','fontsize',12,'interpreter','latex')
ylabel('Power Draw (W)','fontsize',12,'interpreter','latex')
ylim([500,1500])
xlim([0,timeEnd/60])
title('Vision Power Draw per Orbital Period','fontsize',12,'interpreter','latex')
%xline(timeCal/60,'Label','Calibration Mode','LabelHorizontalAlignment','left','LabelVerticalAlignment','top')
%xline(timeNom/60,'Label','Nominal Mode','LabelVerticalAlignment','top')
%xline(timeEcl/60,'Label','Eclipse Mode','LabelVerticalAlignment','top')

