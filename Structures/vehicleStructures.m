% stress analysis on the rover body and orbiter bus
%% Computing the 
maxGs = 8.5; %from falcon guide
maxGravitationalAcceleration = 9.81*maxGs;

% defining rover quantities
betaRover = 0.4974; %area parameter
roverBusLength = .892; %m
roverBusWidth = .462; %m
roverAppliedMass = 65.73; %internal mass of rover, kg
roverAppliedForce = maxGravitationalAcceleration*roverAppliedMass;
roverBusPressureLoad = roverAppliedForce/(roverBusWidth*roverBusLength);

% defining orbiter quantities
betaSpacecraft = .4356; %area parameter
spacecraftBusLength = 1.05; %m
spacecraftBusWidth = .835; %m
spacecraftAppliedMass = 1239; %internal mass of orbiter, kg
spacecraftAppliedForce = maxGravitationalAcceleration*spacecraftAppliedMass;
spacecraftBusPressureLoad = spacecraftAppliedForce/(spacecraftBusLength*spacecraftBusWidth);

aluminumYieldStress = 503e6; %Pa
FoS = 1.25;
reqThickness = @(maxStress, beta1, appliedPressure, width) (sqrt(beta1*appliedPressure*width^2/(maxStress)));
reqThicknessRoverBus = reqThickness(aluminumYieldStress/FoS, betaRover, roverBusPressureLoad, roverBusWidth);
reqThicknessSpacecraftBus = reqThickness(aluminumYieldStress/FoS, betaSpacecraft, spacecraftBusPressureLoad, spacecraftBusWidth);




