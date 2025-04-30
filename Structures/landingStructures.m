% stress analysis on the landing platform and heat shield for EDL
%% Landing platform stresses
% defining constants
aluminumYieldStress = 503e6; %Pa
aluminumDensity = 2810; %kgm/m^3
titaniumYieldStress = 880e6; %Pa
titaniumDensity = 4430; %kg/m^3
roverMass = 425.9; %kg
deltaV = 5; %m/s
impulseTime = 0.1; %s, conservative
length = 1.675; %m
width = 1.25; %m

% area parameter
beta_1 = 0.4356;

% compute applied force/pressure
appliedForce = roverMass*deltaV/impulseTime; %N
appliedPressure = appliedForce/(length*width); %Pa

% func to compute required thickness based on max allowable stress
% from roarks, max stress is along the longitudinal axis, and solving for t
% yields the following
reqThickness = @(maxStress) (sqrt(beta_1*appliedPressure*width^2/(maxStress)));

% use fos of 1.25 to compute maximum allowable stresses, then apply func
FoS = 1.25;
reqThicknessTitanium = reqThickness(titaniumYieldStress/FoS);
reqThicknessAluminum = reqThickness(aluminumYieldStress/FoS);

% compute masses (verified in sw)
massTitanium = reqThicknessAluminum*length*width*titaniumDensity;
massAluminum = reqThicknessTitanium*length*width*aluminumDensity;

%% Internal heat shield structure stresses
% defining constants
entryMass = 1141.77; %kg 
accel = 22; %m/s^2
rb = 1.778; %m
v = 0.34; %poisson's ratio, dimensionless

% computing applied force/pressure
appF = entryMass*accel; %N
area = pi*1.7^2;
pressure = appF/area;

% again using a FoS of 1.25 for allowable stress
FoS = 1.25;
allowableStress = aluminumYieldStress/FoS;
const = 3*(3+v)/8;

% computing thickness and mass
t = sqrt(const*pressure*rb^2/(allowableStress));
mass = pi*rb^2*t*aluminumDensity;