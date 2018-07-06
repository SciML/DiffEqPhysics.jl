abstract type Thermostat
end

struct NullThermostat <: Thermostat
end


struct AndersenThermostat{νType<:Real, tType<:Real, kType<:Real} <: Thermostat
    ν::νType
    T::tType
    kb::kType
end

struct BerendsenThermostat{τType<:Real, tType<:Real} <: Thermostat
    T::tType
    τ::τType
    γ::tType
end

function BerendsenThermostat(T::Real, τ::Real)
    BerendsenThermostat(T, τ, 0.5/τ)
end

function berendsen_acceleration!(dv, v, ms, kb, N, Nc, p::BerendsenThermostat)
    T = md_temperature(v, ms, kb, N, Nc)
    @. dv += p.γ*(p.T/T-1)*v
end

# N - number of particles
# Nc - number of constraints
function md_temperature(vs, ms, kb, N, Nc)
    e_kin = sum(dot(ms, sum(vs.^2,1)))
    temperature = e_kin/(kb*(3*N-Nc))
    return temperature
end