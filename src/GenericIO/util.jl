module util

using PeriodicTable, Unitful

export dalton2elem

function dalton2elem(mass::Real)
    if mass > 239
        error("TRans-Uranium is not supported. Detected: $mass")
    elseif mass <= 0
        error("Atomic mass must be positive. Detected: $mass")
    end
    for e in elements
        if ≈(ustrip(e.atomic_mass), mass, atol=1e-2)
            return Symbol(e.symbol)
        end
    end

    error("Unknown mass of atoms $mass detected. Atomic mass must be in Dalton and be have 1e2 presicion or more.")
end

function dalton2elem(mass::AbstractString)
    mass = parse(Float64, mass)
    dalton2elem(mass)
end

end