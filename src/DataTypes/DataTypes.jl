module DataTypes

using Graphs
using LinearAlgebra
using MetaGraphs
using StaticArrays

import Base: push!, getindex

include("util.jl")
include("HeierchyLabels/HeierchyLabels.jl")

using  .HeierchyLabels

import .HeierchyLabels: add_label!, labels

export Position, BoundingBox, System
export get_position, set_position!
export Id

#####
##### Type `Position` definition
#####

const Position{D, F} = Vector{MVector{D, F}} where {D, F <: AbstractFloat}

function Position(D::Integer, F::Type{<:AbstractFloat}, natom::Integer)
    return if natom < 0
        error("natom must be >= 0")
    elseif natom == 0
        Vector{MVector{D, F}}(undef, 0)
    else
        [MVector{D, F}(zeros(F, D)) for _ in 1:natom]
    end
end

function Position{D, F}() where {D, F}
    Position(D, F, 0)
end

function Position{D, F}(n::Integer) where {D, F}
    Position(D, F, n)
end

function Position()
    Position(3, Float64, 0)
end

function Position(n::Integer)
    Position(3, Float64, n)
end

function push!(p::Position, x::AbstractVector)
    D = length(x)
    F = eltype(x)
    push!(p, MVector{D, F}(x))
end

function set_position!(p::Position, id::Integer, x::AbstractVector)
    p[id] .= x
end
precompile(set_position!, (Position, Int64, Vector{Float64}))

function get_position(p::Position, id::Integer)
    p[id]
end
precompile(get_position, (Position, Int64))

#####
##### Type `BoundingBox` definition
#####

struct BoundingBox{D, F <: AbstractFloat}
    origin::SVector{D, F}
    axis::SMatrix{D, D, F}
    BoundingBox(origin::SVector{D, F}, axis::SMatrix{D, D, F}) where {D, F} = new{D, F}(origin, axis)
end

function BoundingBox(
    D::Integer, F::Type{<:AbstractFloat}, origin::AbstractVector, axis::AbstractMatrix
)
    if D <= 0
        error("D must be positive integer ")
    elseif !(D == length(origin)) == size(axis, 1) == size(axis, 2)
        error("Dimension mismatch ")
    end

    d = det(axis)
    if d < 0
        error("left-handed system not allowed")
    elseif d == 0
        error("box is degenerated")
    end

    return BoundingBox(SVector{D, F}(origin), SMatrix{D, D, F}(axis))
end

function BoundingBox{D, F}(origin::AbstractVector, axis::AbstractMatrix) where {D, F<:AbstractFloat}
    BoundingBox(D, F, origin, axis)
end

function BoundingBox(origin::AbstractVector, axis::AbstractMatrix)
    D = length(origin)
    F = eltype(origin)
    BoundingBox(D, F, origin, axis)
end
precompile(BoundingBox, (Vector{Float64}, Matrix{Float64}))

#####
##### Type `PropMap` definition
#####

mutable struct PropMap
    #atomtype``
end

#####
##### Type `System` definition
#####

mutable struct System{D, F<:AbstractFloat}
    time::F
    topology::Graph{<:Integer}
    box::BoundingBox{D, F}

    # atom property
    position::Position{D, F}
    elem::Vector{String}

    label_heierchy::Dict{Symbol, LabelHeierchy}
    props::Dict{<:AbstractString, PropMap}
end

function select_atom(s::System, label::Label)
    #return {atom id | atom \in Label}
end

function props(s::System, label::Label)
    # return Dict like label => property
end

function atomtype(s::System, id::Integer)
    s.atype[id]
end

function set_atomtype!(s::System, id::Integer, type)
    d.atype[id] = type
end

function add_atomtype!(s::System, type)
    push!(s.atype, type)
end

function elem(s::System, id::Integer)
    s.elem[id]
end

function set_elem!(s::System, id::Integer, elem)
    d.elem[id] = elem
end

function add_elem!(s::System, elem)
    push!(s.elem, elem)
end

function prop(s::System, property::Symbol)
    s.props[property]
end

for field in [:time, :topology, :box, :heierchy, :labels]
    @eval function $field(s::System)
        s.($field)
    end
end

function positions(s::System)
    s.position
end

function get_position(s::System, id::Integer)
    get_position(positions(s), id)
end
precompile(get_position, (System, Int64))

function set_position!(s::System, id::Integer, x::AbstractVector)
    set_position!(positions(s), id, x)
end
precompile(set_position!, (System, Int64, Vector{Float64}))

function add_atom!(s::System; connect::Vector{<:Integer}, pos, type, elem)
    @assert length(positions(s)) == nv(topology(s))
    id = length(positions(s)) + 1

    push!(positions(s), pos)
    add_atomtype!(s, type)
    add_elem!(s, elem)
    for key in keys(labels(s))
        push!(labels(s)[key], label)
    end # add_label!(s, key, value), add_label!(s, :, value)

    t = topology(s)
    id = nv(t) + 1
    if !add_vertex!(t)
        error("atom addition failed at atom id $id")
    end
    for v in connect
        if !add_edge!(t, id, v)
            error("bond addition between $id and $v failed ")
        end
    end
end

function remove_bond(s::System, v1::Integer, v2::Integer)
    println("yet implemented")
end

function set_box!(s, origin::AbstractVector, axis::AbstractMatrix)
    s.box = BoundingBox(origin, axis)
end

function System()
    System{3}()
end
precompile(System, ())

function System{D}() where {D}
    System{D}(
        Vector{Int64}(undef, 0),
        Graph(),
        BoundingBox(),
        Dict{Symbol, LabelHeierchy}(),
        Position{D, Float64}(),
        Vector{Int64}(undef, 0),
        Vector{String}(undef, 0),
        Dict{Symbol, Vector{Int64}}(),
        Dict{Int64, Any}()
    )
end

include("interface.jl")

function System(g::MetaGraph)
    D = props(g, 1)[:position] |> length
    F = props(g, 1)[:position] |> eltype

    s = System()
    s.topology = SimpleGraph(g)
    s.box = BoundingBox()
    #s.elem = map(v -> props(g, v)[:element], vertices(g))

    for v in vertices(g)
        #push!(s.position, MVector{D, F}(props(g, v)[:position]))
        add_atom!(s,
            connect = all_neighbors(g, v),
            pos  = props(g, v)[:position],
            type = 0,
            elem = props(g, v)[:elem])
    end
    s
end

"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""



end #module
