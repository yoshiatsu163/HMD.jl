module DataTypes

export Position, BoundingBox, System
export get_position, set_position!
export Id

using StaticArrays, Graphs, MetaGraphs

import Base: push!, getindex

include("util.jl")
include("HeierchyLabels/HeierchyLabels.jl")
using  .HeierchyLabels
import .HeierchyLabels: add_label!, labels

const Position{D, F} = Vector{MVector{D, F}} where {D, F <: AbstractFloat}
function Position(D, F, natom)
    [MVector{D, F}(zeros(D)) for _ in 1:natom]
end

function Position{D, F}() where {D, F}
    Vector{MVector{D, F}}(undef, 0)
end

function Position()
    Vector{MVector{3, Float64}}(undef, 0)
end

function Base.push!(p::Position, x::AbstractVector)
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

# 関数の引数がパラメタ付き型の場合はwhere {}を使うことで関数内でパラメタアクセスが可能
# e.g. f(a::SVector{D, F}) where {D, F}
struct BoundingBox{D, F <: AbstractFloat}
    origin::SVector{D, F}
    axis::SMatrix{D, D, F}
    #BoundingBox(origin::SVector{D, F}, axis::SMatrix{D, D, F}) where {D, F} = new{D, F}(origin, axis)
end

function BoundingBox(origin::AbstractVector, axis::AbstractMatrix)
    #d = det(axis)
    #if d < 0
    #    error("left-handed system not allowed")
    #elseif d == 0
    #    error("box degenerated")
    #end
    D = length(origin)
    F = eltype(origin)
    BoundingBox(SVector{D, F}(origin), SMatrix{D, D, F}(axis))
end
precompile(BoundingBox, (Vector{Float64}, Matrix{Float64}))

function BoundingBox()
    BoundingBox(zeros(3), zeros(3, 3))
end

mutable struct System{D}
    time::Vector{<:Integer}                       
    topology::Graph{<:Integer}             
    box::BoundingBox                      

    # atom property
    position::Position{D, <:AbstractFloat}        
    atype::Vector{<:Any}                
    elem::Vector{String}

    heierchy::Dict{Symbol, LabelHeierchy}
    props::Dict{<:Integer, <:Any}
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

function add_atom!(s::System; connect::Vector{<:Integer} pos, type, elem)
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