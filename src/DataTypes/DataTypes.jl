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

import Base: getindex
import Base: >, <, >=, <=, +, -, *, /, ==, string

include("util.jl")
include("HierarchyLabels/HierarchyLabels.jl")

using  .HierarchyLabels

export Position, BoundingBox, AbstractSystem, System, Label
export time, set_time!, natom, topology, box, set_box!
export all_element, element, set_element!
export all_positions, position, set_position!
export hierarchy_names, hierarchy, add_hierarchy!, remove_hierarchy!, merge_hierarchy!
export prop_names, props, prop, labels_in_prop, add_prop!, set_prop!
export labels, add_label!, add_relation!, insert_relation!, remove_label!, remove_relation!
export Id, Category

export contains, has_relation, issuper, issub, super, sub

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


#####
##### Type `System` definition
#####

abstract type AbstractSystem end

mutable struct System{D, F<:AbstractFloat} <: AbstractSystem
    time::F
    topology::Graph{<:Integer}
    box::BoundingBox{D, F}

    # atom property
    position::Position{D, F}
    element::Vector{Category{Element}}

    hierarchy::Dict{<:AbstractString, LabelHierarchy}
    props::Dict{<:AbstractString, Dict{Label, Any}}
end

include("property.jl")

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

function natom(s::AbstractSystem)
    length(s.position)
end

function time(s::AbstractSystem)
    s.time
end

function set_time!(s::AbstractSystem, time::AbstractFloat)
    s.time = time
end

function topology(s::AbstractSystem)
    s.topology
end

function box(s::AbstractSystem)
    s.box
end

function set_box!(s::AbstractSystem, box::BoundingBox)
    s.box = box
end

function all_element(s::AbstractSystem)
    s.element
end

function element(s::AbstractSystem, atom_id::Integer)
    s.element[atom_id]
end

function set_element!(s::AbstractSystem, atom_id::Integer, ename::AbstractString)
    s.element[atom_id] = Category{Element}(ename)
end

function set_element!(s::AbstractSystem, atom_ids::AbstractVector{<:Integer}, enames::AbstractVector{<:AbstractString})
    if length(atom_ids) != length(enames)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.element[atom_ids] .= Category{Element}.(enames)
end

function all_positions(s::AbstractSystem)
    s.position
end

function position(s::AbstractSystem, atom_id::Integer)
    s.position[atom_id]
end

function set_position!(s::AbstractSystem, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] .= x
end

function set_position!(s::AbstractSystem, atom_ids::AbstractVector{<:Integer}, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    if length(atom_ids) != length(x)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.position[atom_ids] .= x
end

function hierarchy_names(s::AbstractSystem)
    s.hierarchy |> keys
end

function hierarchy(s::AbstractSystem, hname::AbstractString)
    s.hierarchy[hname]
end

function merge_hierarchy!(s::AbstractSystem, hie1::Dict{<:AbstractString, LabelHierarchy}, hie2::Dict{<:AbstractString, LabelHierarchy})
    s.hierarchy = merge(hie1, hie2)
end

function add_hierarchy!(s::AbstractSystem, hname::AbstractString)
    if hname in hierarchy_names(s)
        error("hierarchy $(hname) already exists. ")
    end
    push!(s.hierarchy, hname => LabelHierarchy())
    _add_label!(hierarchy(s, hname), Entire_System)
end

function remove_hierarchy!(s::AbstractSystem, hname::AbstractString)
    delete!(s.hierarchy, hname)
end

function labels(s::AbstractSystem, hname::AbstractString)
    lh = hierarchy(s, hname)
    return _labels(lh)
end

function add_label!(s::AbstractSystem, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)

    result = _add_label!(lh, label)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return nothing
        _              => error("fatal error")
    end
end

function add_relation!(s::AbstractSystem, hname::AbstractString; super::Label, sub::Label)
    lh = hierarchy(s, hname)

    result = _add_relation!(lh; super=super, sub=sub)
    @match result begin
        Label_Missing     => error("Super or sub not found. ")
        Label_Duplication => error("super and sub are equal. ")
        Relation_Occupied => error("""There in already relation between super and sub labels. Please use "insert_relation!()" instead. """)
        success           => return nothing
        _                 => error("fatal error")
    end
end

function insert_relation!(s::AbstractSystem, hname::AbstractString, label::Label; super::Label, sub::Label)
    lh = hierarchy(s, hname)

    add_label!(s, hname, label)
    add_relation!(s, hname; super=super, sub=label)
    add_relation!(s, hname; super=label, sub=sub)
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function remove_label!(s::AbstractSystem, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)

    if !_contains!(lh, label)
        error("label $(label) not found. ")
    end
    @assert _remove_label!(lh, label)

    return nothing
end

function remove_relation!(s::AbstractSystem, hname::AbstractString; super::Label, sub::Label)
    lh = hierarchy(s, hname)

    if !_has_relation(lh, super, sub)
        error("relation betewwn $(super) and $(sub) not found. ")
    end
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function contains(s::AbstractSystem, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)
    return _contains(lh, label)
end

function issuper(s::AbstractSystem, hname::AbstractString, label1::Label, label2::Label)
    lh = hierarchy(s, hname)
    return _issuper(lh, label1, label2)
end

function issub(s::AbstractSystem, hname::AbstractString, label1::Label, label2::Label)
    lh = hierarchy(s, hname)
    return _issub(lh, label1, label2)
end

function super(s::AbstractSystem, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)
    return _super(lh, label)
end

function sub(s::AbstractSystem, hname::AbstractString, label::Label)
    lh = hierarchy(s, hname)
    return _sub(lh, label)
end

include("test.jl")

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
