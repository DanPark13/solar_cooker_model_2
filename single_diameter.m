function [elapsedTime, waterTemperature] = single_diameter()
    %% DESCRIPTION sweep_diameter 
    % Performs a parameter sweep over the diameter of the mirror of the parabolic solar cooker
    
    %% TODO: add pot temperature change, add convection, add radiation, add in diameter sweep
    
    %% Variables
    initialTemperature = 293.15;    % K
    airTemperature = 273.15;        % K (Himalayas) % to be simulated
    waterVolume = 2;                % L
    
    % Radiation
    boltzmannConstant = 5.67 * 10^-8;    % W/m^2K^4
    efficiency_of_Absorption = 1; % Unitless
    insolation = 1000; % W/m^2
    mirror_diameter = 1; % m
    mirror_area = pi * ((mirror_diameter/2)^2);
     
    % pot dimensions refer to inner dimensions
    potDiameter = 14/100;           % m
    potThickness = .5/100;          % m
    potHeight = 14/100;             % m
    potConductivity = 5;          % thermal conductivity of pot walls (W/(m*K))
    potSpecificHeat = 380;          % specific heat of water (J/kg*K)
    waterSpecificHeat = 4186;       % specific heat of water (J/kg*K)
    waterDensity = 1;               % kg/L
    
    %% Calculated Variables
    potOuterDiameter = potDiameter + 2*potThickness;
    potOuterHeight = potHeight + 2*potThickness;
    
    potRadius = potDiameter/2;
    waterHeight = .001 * waterVolume / (pi*potRadius^2);
    conductionArea = pi * potDiameter * waterHeight + potRadius^2;
    
    convectionArea = 2*pi*(potDiameter/2)^2 + pi*potDiameter*potOuterHeight;
    waterMass = waterDensity * waterVolume;
    initialEnergy = temperatureToEnergy(initialTemperature, waterMass, waterSpecificHeat);
    
    %% Simulation Parameters
    timeParams = [0, 30*60];   % convert the minutes to seconds
    
    %% Simulation 
    [T, U] = ode45(@heatFlow, timeParams, initialEnergy);

    elapsedTime = T./60;
    waterTemperature = energyToTemperature(U, waterMass, waterSpecificHeat);
    % potTemperature = energyToTemperature(U, potMass, potSpecificHeat);
    % environmentLost
    
    function dUdt = heatFlow(~, U)
        deltaTemp = energyToTemperature(U, waterMass, waterSpecificHeat) - airTemperature;
        radiationRate = efficiency_of_Absorption * insolation * mirror_area;
        % conductionRate = potConductivity * conductionArea / potThickness * deltaTemp;
        convectionRate = 25*convectionArea * deltaTemp;
        dUdt =  radiationRate - convectionRate;%conductionRate;
    end
end