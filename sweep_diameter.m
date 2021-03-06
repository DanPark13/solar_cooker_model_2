function [TimeToBoil,Diameter] = sweep_diameter(diameterRange, windSpeed, airTemperature, sunlightIntensity)
    %% DESCRIPTION sweep_diameter 
    % Performs a parameter sweep over the diameter of the mirror of the parabolic solar cooker
    
    %% TODO: add pot temperature change, add convection, add radiation, add in diameter sweep
    
    %% Variables
    initialTemperature = 293.15;    % K
    waterVolume = 2;                % L
    
    % pot dimensions refer to inner dimensions
    potDiameter = 14/100;           % m
    potThickness = .5/100;          % m
    potHeight = 14/100;             % m
    potConductivity = 401;          % thermal conductivity of pot walls (W/(m*K))
    potSpecificHeat = 380;          % specific heat of water (J/kg*K)
    waterSpecificHeat = 4186;       % specific heat of water (J/kg*K)
    waterDensity = 1000;            % kg/m^3
    
    %% Calculated Variables
    potOuterDiameter = potDiameter + 2*potThickness;
    potOuterHeight = potHeight + 2*potThickness;
    
    potRadius = potDiameter/2;
    waterHeight = waterVolume / (pi*potRadius^2);
    conductionArea = pi * potDiameter * waterHeight + potRadius^2;
    
    convectionArea = 2*pi*(potOuterDiameter/2)^2 + pi*potOuterDiameter*potOuterHeight; 
    waterMass = waterDensity * waterVolume;
    initialEnergy = temperatureToEnergy(initialTemperature, waterMass, waterSpecificHeat);
    
    %% Simulation Parameters
    timeParams = [0, 30 * 60];   % convert the minutes to seconds
    
    %% Simulation 
    [T, U] = ode45(@heatFlow, timeParams, initialEnergy);
    
    elapsedTime = T./60;
    waterTemperature = energyToTemperature(U, waterMass, waterSpecificHeat);
    potTemperature = energyToTemperature(U, potMass, potSpecificHeat);
    % environmentLost
    
    % Radiation Plotting
    plot(elapsed_time, potTemperature, '*-'); % Add Convection Variables
    xlabel('Time (minutes)');
    ylabel('Temperature (K)');
    title('Convection Over Time')
    
    % Conduction Plotting
    plot(elapsed_time, waterTemperature, '*-');
    xlabel('Time (minutes)');
    ylabel('Temperature (K)');
    title('Water Boiling Over Time')
    
    % Convection Plotting
    % TO_DO
    
    function dUdt = heatFlow(~, U)
        deltaTemp = energyToTemperature(U, waterMass, waterSpecificHeat) - airTemperature;
        dUdt = -cupConductivity * conductionArea / cupThickness * deltaTemp;
    end
    
end