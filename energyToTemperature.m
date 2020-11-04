%% res = energyToTemperature(U, m, c)
% Convert internal energy (J) to temperature (K) for a given mass (kg)
% of a substance with a specific heat (J/(kg*K)).
function res = energyToTemperature(U, m, c)
    res = U / heatCapacity(m, c);
end