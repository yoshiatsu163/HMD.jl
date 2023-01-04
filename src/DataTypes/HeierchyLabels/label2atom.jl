mutable struct Label2atom{T1, T2}
    toatom::Dict{T1, T2}
    tolabel::Vector{Vector{Label}}
end

function Label2atom()
    Label2atom{Label, Vector{Int64}}(
        Dict{Label, Vector{Int64}}(), Vector{Vector{Label}}())
end

function toatom(l2a::Label2atom)
    l2a.toatom
end

function tolabel(l2a::Label2atom)
    l2a.tolabel
end

function _labels(l2a::Label2atom)
    toatom(l2a) |> keys |> collect
end

function getindex(l2a::Label2atom, label::Label)
    toatom(l2a)[label]
end

function getindex(l2a::Label2atom, atomid::Integer)
    tolabel(l2a)[atomid]
end

# add labelのみだとlabelのついてない原子を扱えない
# 1. Heierchyにそもそも入れない
# 2. add_atom等を用意してatomとのコヒーレンシを保つ
#   sync!()でもいいがテストユニットが大きくなる
function _add_label!(l2a::Label2atom, label::Label, ids::Vector{<:Integer})
    if label ∈ toatom(l2a) |> keys
        error("label already exists in Label => atom mapping. ")
    end

    push!(toatom(l2a), label => ids)
    for id in ids
        push!(tolabel(l2a)[id], label)
    end
end

function _add_atom!(l2a::Label2atom, n::Integer)
    append!(tolabel(l2a), [Label[] for _ in 1:n])
    tolabel(l2a) |> length
end
