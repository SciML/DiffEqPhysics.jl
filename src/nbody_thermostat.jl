abstract type Thermostat
end

struct NullThermostat <: Thermostat
end


struct AndersenThermostat{νType<:Real, tType<:Real, kType<:Real} <: Thermostat
    ν::νType
    T::tType
    kb::kType
end