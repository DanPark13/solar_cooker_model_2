%% res = heatCapacity(mass, specificHeat)
% Compute the heat capacity (J/K) of a given mass (kg) of a substance
% with a given specific heat (J/(kg*K)).
function res = heatCapacity(mass, specificHeat)
    res = mass * specificHeat;
end