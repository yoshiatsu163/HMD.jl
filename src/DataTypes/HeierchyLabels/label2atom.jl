mutable struct Label2atom{T1, T2}
    forward::Dict{T1, T2}
    reverse::Vector{Vector{Label}}
end

function Label2atom()
    Label2atom{Label, Int64}(Dict{Label, Int64}(), Label[][])
end

function forward(l2a::Label2atom)
    l2a.forward
end

function revesre(l2a::Label2atom)
    l2a.reverse
end

function _labels(l2a::Label2atom)
    forward(l2a) |> keys
end

function getindex(l2a::Label2atom, label::Label)
    forward(l2a)[label]
end

function getindex(l2a::Label2atom, atomid::Integer)
    reverse(l2a)[atomid]
end

function _add_label!(l2a::Label2atom, label::Label, ids::Vector{<:Integer})
    push!(forward(l2a), label => ids)
    for id in ids
        push!(l2a[id], label)
    end
end
