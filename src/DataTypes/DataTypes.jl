# presetとしてpolymer hierarchyを追加: polymer >: {monomer >: connection, endcap >: endatom}
# add_label, add, \oplus, add!を作成 テスト可能
# selectionなど読み込み機能を追加
#
#

module DataTypes

using Graphs
using LinearAlgebra
using Match
using MetaGraphs
using PeriodicTable
using StaticArrays

import Base: push!, getindex
import Base: >, <, >=, <=, +, -, *, /, ==, string

include("util.jl")
include("HierarchyLabels/HierarchyLabels.jl")

using  .HierarchyLabels

export Position, BoundingBox, System
export get_position, set_position!
export Id, Category

const Entire_System = Label(1, "entire_system")

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

function BoundingBox{D, F}() where {D, F<:AbstractFloat}
    BoundingBox{D, F}(zeros(F, 3), Matrix{F}(I, D, D))
end

function BoundingBox(origin::AbstractVector, axis::AbstractMatrix)
    D = length(origin)
    F = eltype(origin)
    BoundingBox(D, F, origin, axis)
end
precompile(BoundingBox, (Vector{Float64}, Matrix{Float64}))

#####
##### Type `System` definition
#####

mutable struct System{D, F<:AbstractFloat}
    time::F
    topology::Graph{<:Integer}
    box::BoundingBox{D, F}

    # atom property
    position::Position{D, F}
    element::Vector{Category{Element}}

    hierarchy::Dict{<:AbstractString, LabelHierarchy}
    props::Dict{<:AbstractString, Dict{Label, Any}}
end

function System{D, F}() where {D, F<:AbstractFloat}
    System{D, F}(
        zero(F),
        SimpleGraph{Int64}(),
        BoundingBox{D, F}(),
        Position{D, F}(),
        Vector{Category{Element}}(undef, 0),
        Dict{String, LabelHierarchy}(),
        Dict{String, Dict{Label, Any}}()
    )
end

function natom(s::System)
    length(s.position)
end

function time(s::System)
    s.time
end

function set_time!(s::System, time::AbstractFloat)
    s.time = time
end

function topology(s::System)
    s.topology
end

function box(s::System)
    s.box
end

function set_box!(s::System, box::BoundingBox)
    s.box = box
end

function all_element(s::System)
    s.element
end

function element(s::System, atom_id::Integer)
    s.element[atom_id]
end

function set_element!(s::System, atom_id::Integer, ename::AbstractString)
    s.element[atom_id] = Category{Element}(ename)
end

function set_element!(s::System, atom_ids::AbstractVector{<:Integer}, enames::AbstractVector{<:AbstractString})
    if length(atom_ids) != length(enames)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.element[atom_ids] .= Category{Element}.(enames)
end

function all_positions(s::System)
    s.position
end

function position(s::System, atom_id::Integer)
    s.position[atom_id]
end

function set_position!(s::System, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] .= x
end

function set_position!(s::System, atom_ids::AbstractVector{<:Integer}, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    if length(atom_ids) != length(x)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.position[atom_ids] .= x
end

function hierarchy_names(s::System)
    s.hierarchy |> keys
end

function hierarchy(s::System, hname::AbstractString)
    s.hierarchy[hname]
end

function add_hierarchy!(s::System, hname::AbstractString)
    if hname in hierarchy_names(s)
        error("hierarchy $(hname) already exists. ")
    end
    push!(s.hierarchy, hname => LabelHierarchy())
    _add_label!(hierarchy(s, hname), Entire_System, super=No_Label, sub=No_Label)
end

function labels(s::System, hname::AbstractString)
    lh = hierarchy(s, hname)
    return _labels(lh)
end

function add_label!(s::System, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)

    if label ∈ lh
        error("Label $(label) already exists. ")
    end

    _add_label!(lh, label, super=No_Label, sub=No_Label)

    return nothing
end

function add_relation!(s::System, hname::AbstractString; super::Label, sub::Label)
    lh = hierarchy(s, hname)

    if super ∉ lh || sub ∉ lh
        error("Super or sub not found. ")
    elseif _has_relation(lh, super, sub)
        error("""There in already relation between super and sub labels. Please use "insert_relation!()" instead. """)
    end

    return nothing
end

function insert_relation!(s::System, hname::AbstractString, label::Label; super::Label, sub::Label)
    lh = hierarchy(s, hname)

    if label ∈ lh
        error("Label $(label) already exists. ")
    elseif super ∉ lh || sub ∉ lh
        error("Super or sub not found. ")
    elseif !_has_relation(lh, super, sub)
        error("""There is no relation between super and sub. """)
    end

    add_label!(s, hname, label)
    @assert _remove_relation!(lh, super, sub)
    add_relation!(s, hname; super=super, sub=label)
    add_relation!(s, hname; super=label, sub=sub)

    return nothing
end

function prop_names(s::System)
    s.props |> keys
end

function props(s::System, pname::AbstractString)
    s.props[pname]
end

function prop(s::System, pname::AbstractString, label::Label)
    prop(s, pname)[label]
end

function labels_in_prop(s::System, pname::AbstractString)
    props(s, pname) |> keys
end

function ∈(label::Label, prp::Dict{Label, Any})
    label ∈ keys(prp)
end

function ∋(prp::Dict{Label, Any}, label::Label)
    label ∈ keys(prp)
end

function ∉(label::Label, prp::Dict{Label, Any})
    label ∉ keys(prp)
end

function add_prop!(s::System, pname::AbstractString)
    push!(s.props, pname => Dict{Label, Any}())
end

function add_prop!(s::System, pname::AbstractString, label::Label, p::Any)
    prp = props(s, pname)
    if label ∈ prp
        error("label $(label) already exists. ")
    end
    push!(prp, label => p)
end

function set_prop!(s::System, pname::AbstractString, label::Label, p::Any)
    if label ∉ props(s, pname)
        error("label $(label) not found in property $(pname). ")
    end
    props(s, pname)[label] = p
end

function add_atom!(s::System, x::AbstractVector{<:AbstractFloat}, elem::AbstractString)
    atom_id = natom(s) + 1
    push!(s.position, x)
    push!(s.element , Category{Element}(elem))
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        lh = hierarchy(s, hname)
        _add_label!(lh, atom_label(atom_id), super=Entire_System, sub=No_Label)
    end

    return nothing
end

function add_bond!(s::System, atom_id1::Integer, atom_id2::Integer)
    topo = topology(s)
    @assert add_edge!(topo, atom_id1, atom_id2)

    return nothing
end

function atom_label(atom_id::Integer)
    Label(atom_id, "")
end



include("interface.jl")

#function System(g::MetaGraph)
#    D = props(g, 1)[:position] |> length
#    F = props(g, 1)[:position] |> eltype
#
#    s = System()
#    s.topology = SimpleGraph(g)
#    s.box = BoundingBox()
#    #s.elem = map(v -> props(g, v)[:element], vertices(g))
#
#    for v in vertices(g)
#        #push!(s.position, MVector{D, F}(props(g, v)[:position]))
#        add_atom!(s,
#            connect = all_neighbors(g, v),
#            pos  = props(g, v)[:position],
#            type = 0,
#            elem = props(g, v)[:elem])
#    end
#    s
#end

"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""



end #module
