%% res = temperatureToEnergy(T, m, c)
% Convert temperature (K) to internal energy (J) for a given mass (kg)
% of a substance with a specific heat (J/(kg*K)).
function res = temperatureToEnergy(T, m, c)
    res = T * heatCapacity(m, c);
end