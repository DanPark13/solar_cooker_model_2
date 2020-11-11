function [T, waterTemperature, potTemperature, endTime] = single_diameter_v2(potMaterial, mirrorDiameter)
    %% sweep_diameter_v2
    %  Simulates 2 liters of water in a pot suspended in a parabolic solar
    %  cooker. The simulation ends when the water boils, or when an hour
    %  has passed.
    %
    %  returns:
    %   T: column vector of each time step from the simulation
    %   waterTemperature: column vector of the water temperature at each
    %       time step
    %   potTemperature: column vector of the pot's temperature at each
    %       time step
    %   endTime: the number of seconds passed before the water boils
    
    %% Variables
    waterVolume = 2;    % L
    
    % material properties
    copper = struct(...
        'density', 8920,...         % kg/m^3
        'specificHeat', 380,...     % J/kg K
        'conductivity', 401);       % W/m K
    
    aluminum = struct(...
        'density', 2700,...         % kg/m^3
        'specificHeat', 900,...     % J/kg K
        'conductivity', 237);       % W/m K
    
    stainlessSteel304 = struct(...
        'density', 8030,...         % kg/m^3
        'specificHeat', 500,...     % J/kg K
        'conductivity', 16.2);      % W/m K
    
    water = struct(...
        'density', 1,...            % kg/L
        'specificHeat', 4186);      % J/kg K
    
    % pot properties
    switch potMaterial
        case 'Copper'
            potMaterial = copper;
        case '304 Stainless Steel'
            potMaterial = stainlessSteel304;
        case 'Aluminum'
            potMaterial = aluminum;
        otherwise
            potMaterial = stainlessSteel304;
    end
    potInnerDiameter = 14/100;      % m
    potThickness = .5/100;          % m
    potDepth = 14/100;              % m
    efficiencyOfAbsorption = .75;   % Unitless
    
    % environment
    initialTemperature = 273.15;    % K
    airTemperature = 243.15;        % K (Himalayas) % to be simulated
    insolation = 921;              % W/m^2 (solar noon at 28 degrees latitude on January 1)
    convectionCoefficient = 25;     % W/m^2*K, could range from 2-25, depending on wind speed
    

    %% Calculated Variables
    % pot properties
    potOuterDiameter = potInnerDiameter + 2*potThickness;   % m
    potOuterHeight = potDepth + 2*potThickness;       % m
    potRadius = potInnerDiameter/2; % m
    potVolume = pi*(potOuterDiameter/2)^2*potOuterHeight - pi*potRadius^2*potDepth;   % includes lid
    potMass = potVolume*potMaterial.density;
    
    % mirror properties
    mirror_area = pi * ((mirrorDiameter/2)^2); % m^2
    
    % energy transfer calculations
    waterHeight = .001 * waterVolume / (pi*potRadius^2);    % m, .001 is to convert L to m^3
    waterMass = water.density * waterVolume;                % kg
    
    conductionArea = pi * potInnerDiameter * waterHeight + potRadius^2; % m^2
    convectionArea = 2*pi*(potOuterDiameter/2)^2 + pi*potOuterDiameter*potOuterHeight;  % m^2
    
    initialEnergyWater = temperatureToEnergy(initialTemperature, waterMass, water.specificHeat); % J
    initialEnergyPot = temperatureToEnergy(initialTemperature, potMass, potMaterial.specificHeat);    % J
    
    %% Simulation
    timeParams = 0 : .5 : 24*60*60;   % .5 second time steps, simulation timeout is 1 day
    options = odeset('Events', @eventFunc); % event handle to cut off ode when boiling
    
    [T, U, te, ~, ~] = ode45(@heatFlow, timeParams, [initialEnergyWater, initialEnergyPot], options);
    
    % convert times to minutes before returning
    endTime = te/60;
    T = T./60;
    
    % convert energies to temperatures before returning
    waterTemperature = energyToTemperature(U(:,1), waterMass, water.specificHeat);
    potTemperature = energyToTemperature(U(:,2), potMass, potMaterial.specificHeat);
    
    function [value, isterminal, direction] = eventFunc(~, U)
        value = (U(1) >= temperatureToEnergy(373.15, waterMass,...
            water.specificHeat));    % cuts off ode45 when temp is greater than boiling (373.15 K)
        isterminal = 1;
        direction = 0;
    end
    
    function res = heatFlow(~, U)
        waterTemp = energyToTemperature(U(1), waterMass, water.specificHeat);
        potTemp = energyToTemperature(U(2), potMass, potMaterial.specificHeat);
        
        deltaTempConduction = potTemp - waterTemp;
        conductionRate = potMaterial.conductivity * conductionArea / potThickness * deltaTempConduction;
        
        deltaTempConvection = potTemp - airTemperature;
        convectionRate = convectionCoefficient*convectionArea * deltaTempConvection;
        
        radiationRate = efficiencyOfAbsorption * insolation * mirror_area;
        
        waterRate = conductionRate;
        potRate = radiationRate - convectionRate;
        res = [waterRate; potRate];
    end
end