using DataStructures

Base.@kwdef mutable struct Bimap{D <: AbstractDict, X <: Any, I <: Integer}
    value::Vector{X} = Vector{X}(undef, 0)
    revmap::D{X, I} = Dict()
end

b = Bimap()
