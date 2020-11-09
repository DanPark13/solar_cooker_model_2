function [elapsedTime, waterTemperature, potTemperature] = single_diameter_v2()
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
    copperDensity = .892;            % kg/m^3
    potDiameter = 14/100;           % m
    potThickness = .5/100;          % m
    potHeight = 14/100;             % m
    potConductivity = 5;            % thermal conductivity of pot walls (W/(m*K))
    potSpecificHeat = 380;          % specific heat of water (J/kg*K)
    waterSpecificHeat = 4186;       % specific heat of water (J/kg*K)
    waterDensity = 1;               % kg/L
    
    %% Calculated Variables
    potOuterDiameter = potDiameter + 2*potThickness;
    potOuterHeight = potHeight + 2*potThickness;
    
    potRadius = potDiameter/2;
    
    potVolume = pi*(potOuterDiameter/2)^2*potOuterHeight - pi*potRadius^2*potHeight;
    potMass = potVolume*copperDensity;
    
    waterHeight = .001 * waterVolume / (pi*potRadius^2);
    conductionArea = pi * potDiameter * waterHeight + potRadius^2;
    
    convectionArea = 2*pi*(potOuterDiameter/2)^2 + pi*potOuterDiameter*potOuterHeight;
    waterMass = waterDensity * waterVolume;
    initialEnergyWater = temperatureToEnergy(initialTemperature, waterMass, waterSpecificHeat);
    initialEnergyPot = temperatureToEnergy(initialTemperature, potMass, potSpecificHeat);
    
    %% Simulation Parameters
    timeParams = [0, 10*60];   % convert the minutes to seconds
    initialValues = [initialEnergyWater, initialEnergyPot];
    
    %% Simulation 
    [T, U] = ode45(@heatFlow, timeParams, initialValues);

    elapsedTime = T./60;
    waterTemperature = energyToTemperature(U(:,1), waterMass, waterSpecificHeat);
    potTemperature = energyToTemperature(U(:,2), potMass, potSpecificHeat);
    
    function res = heatFlow(~, U)
        potTemp = energyToTemperature(U(2), potMass, potSpecificHeat);
        waterTemp = energyToTemperature(U(1), waterMass, waterSpecificHeat);
        deltaTempConvection = potTemp - airTemperature;
        deltaTempConduction = potTemp - waterTemp;
        radiationRate = efficiency_of_Absorption * insolation * mirror_area;
        conductionRate = potConductivity * conductionArea / potThickness * deltaTempConduction;
        convectionRate = 25*convectionArea * deltaTempConvection;
        waterRate = conductionRate;
        potRate = radiationRate - convectionRate;
        res = [waterRate; potRate];
    end
end