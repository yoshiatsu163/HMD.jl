# presetとしてpolymer hierarchyを追加: polymer >: {monomer >: connection, endcap >: endatom}
# add_label, add, \oplus, add!を作成 テスト可能
# selectionなど読み込み機能を追加
#
#

module DataTypes

using DataStructures
using Graphs
using LinearAlgebra
using MetaGraphs
using MLStyle
using PeriodicTable
using SimpleWeightedGraphs
using StaticArrays

import Base: getindex, setproperty!, iterate, length
import Base: >, <, >=, <=, +, -, *, /, ==, string, show, convert
import Base: position, time, contains, show
import Base: promote_rule, promote_type
import Base: ∈, ∉
import MetaGraphs: set_prop!, props

include("util.jl")
include("HierarchyLabels/HierarchyLabels.jl")

using  .HierarchyLabels

# core subtype signature
export Position, BoundingBox, HLabel, LabelHierarchy, Id, Category
export >, <, >=, <=, +, -, *, /, ==, string, show, convert, getindex, convert
export id, type, ==, promote_rule, promote_type, length

# core immut signature
export AbstractSystemType, GeneralSystem, AbstractSystem, System
export contains
export all_elements, element
export time, natom, nbond, topology, box, dimension, show
export all_positions, position, travel, wrapped
export hierarchy_names, hierarchy
export all_labels, count_label, has_relation, issuper, issub, super, sub, print_to_string
export prop_names, props, prop, labels_in_prop

# core mut signature
export set_time!, set_box!
export _add_element!, set_element!, _add_position!, set_position!, set_travel!, _change_wrap!
export add_hierarchy!, remove_hierarchy!
export add_prop!, set_prop!
export add_label!, add_labels!, add_relation!, insert_relation!, add_relations!, remove_label!, remove_relation!

# trajectory specific signature

#constants
export Entire_System

const Entire_System = HLabel("entire_system", 1)

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

function Position{D, F}(x::Matrix{F}) where {D, F}
    if size(x, 1) != D
        error("Dimension mismatch")
    end
    return [MVector{D, F}(x[:, i]) for i in 1:size(x, 2)]
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

abstract type AbstractSystem{D, F<:AbstractFloat} end
abstract type AbstractSystemType end

struct GeneralSystem <: AbstractSystemType end

mutable struct System{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractSystem{D, F}
    time::F
    topology::SimpleWeightedGraph{<:Integer, <:Rational}
    box::BoundingBox{D, F}

    # atom property
    position::Position{D, F}
    travel::Vector{MVector{D, Int16}}
    wrapped::Bool
    element::Vector{Category{Element}}

    hierarchy::Dict{<:AbstractString, LabelHierarchy}
    props::Dict{<:AbstractString, Dict{HLabel, Any}}
end

include("property.jl")
include("trajectory.jl")

function System{D, F, SysType}() where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    System{D, F, SysType}(
        zero(F),
        SimpleWeightedGraph{Int64, Rational{Int8}}(),
        BoundingBox{D, F}(),
        Position{D, F}(),
        Vector{MVector{D, Int16}}(undef, 0),
        false,
        Vector{Category{Element}}(undef, 0),
        Dict{String, LabelHierarchy}(),
        Dict{String, Dict{HLabel, Any}}()
    )
end

function System{D, F}() where {D, F<:AbstractFloat}
    return System{D, F, GeneralSystem}()
end

function dimension(s::AbstractSystem{D, F}) where {D, F<:AbstractFloat}
    return D
end

function Base.show(io::IO, ::MIME"text/plain", s::AbstractSystem{D, F}) where {D, F}
    "System{$D, $F}
        time: $(time(s))
        bbox: $(box(s))
        natoms: $(natom(s))
        hierarchy: $(hierarchy_names(s))
    " |> println
end

function Base.show(io::IO, ::MIME"text/plain", s::System{D, F, SysType}) where {D, F, SysType}
    "System{$D, $F, $SysType}
        time: $(time(s))
        bbox: $(box(s))
        natoms: $(natom(s))
        hierarchy: $(hierarchy_names(s))
    " |> println
end

function natom(s::AbstractSystem)
    natm = length(s.position)
    @assert natm == length(s.element)
    return length(s.position)
end

function nbond(s::AbstractSystem)
    return topology(s) |> ne
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

function all_elements(s::AbstractSystem)
    s.element
end

function element(s::AbstractSystem, atom_id::Integer)
    s.element[atom_id]
end

function element(s::AbstractSystem, label::HLabel)
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    s.element[atom_id]

    return nothing
end

function _add_element!(s::AbstractSystem, ename::AbstractString)
    _add_element!(s, Category{Element}(ename))
end

function _add_element!(s::AbstractSystem, ename::Category{Element})
    push!(s.element, ename)
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

function _add_position!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat})
    push!(s.position, x)
    if wrapped(s)
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    push!(s.travel, zeros(Int16, 3))
end

function set_position!(s::AbstractSystem, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] .= x
end

function set_position!(s::AbstractSystem, label::HLabel, x::AbstractVector{<:AbstractFloat})
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    set_position!(s, label, x)

    return nothing
end

function travel(s::AbstractSystem, atom_id::Integer)
    return s.travel[atom_id]
end

function set_travel!(s::AbstractSystem, atom_id::Integer, n::AbstractVector{<:Integer})
    s.travel[atom_id] .= n
end

function wrapped(s::AbstractSystem)
    s.wrapped
end

function _change_wrap!(s::AbstractSystem)
    s.wrapped = !(s.wrapped)
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

function add_hierarchy!(s::AbstractSystem, hname::AbstractString)
    if hname in hierarchy_names(s)
        error("hierarchy $(hname) already exists. ")
    end
    push!(s.hierarchy, hname => LabelHierarchy())
    #_add_label!(hierarchy(s, hname), Entire_System)
    add_label!(s, hname, Entire_System)

    return nothing
end

function remove_hierarchy!(s::AbstractSystem, hname::AbstractString)
    delete!(s.hierarchy, hname)
end

function all_labels(s::AbstractSystem, hname::AbstractString)
    lh = hierarchy(s, hname)
    return _labels(lh)
end

function add_label!(s::AbstractSystem, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    result = _add_label!(lh, label)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return nothing
        _              => error("fatal error")
    end
end

function add_label!(s::AbstractSystem, hname::AbstractString, label_type::Category{HLabel})
    lh = hierarchy(s, hname)

    n = count_label(s, hname, label_type)
    addend = HLabel(label_type, n+1)
    result = _add_label!(lh, addend)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return addend
        _              => error("fatal error")
    end
end

function add_labels!(s::AbstractSystem, hname::AbstractString, label_types::AbstractVector{<:HLabel})
    lh = hierarchy(s, hname)

    _add_labels!(lh, label_types)

    return nothing
end

function count_label(s::AbstractSystem, hname::AbstractString, label_type::Category{HLabel})
    lh = hierarchy(s, hname)
    labels = _label2node(lh) |> keys
    return count(l -> type(l)==label_type, labels)
end

function add_relation!(s::AbstractSystem, hname::AbstractString; super::HLabel, sub::HLabel)
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

function add_relations!(s::AbstractSystem, hname::AbstractString; super::HLabel, subs::AbstractVector{HLabel})
    lh = hierarchy(s, hname)
    for sub in subs
        _add_relation!(lh; super=super, sub=sub, unsafe=true)
    end

    return nothing
end

function insert_relation!(s::AbstractSystem, hname::AbstractString, label::HLabel; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    add_label!(s, hname, label)
    add_relation!(s, hname; super=super, sub=label)
    add_relation!(s, hname; super=label, sub=sub)
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function remove_label!(s::AbstractSystem, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    if !_contains!(lh, label)
        error("label $(label) not found. ")
    end
    @assert _remove_label!(lh, label)

    return nothing
end

function remove_relation!(s::AbstractSystem, hname::AbstractString; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    if !_has_relation(lh, super, sub)
        error("relation betewwn $(super) and $(sub) not found. ")
    end
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function contains(s::AbstractSystem, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _contains(lh, label)
end

function issuper(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issuper(lh, label1, label2)
end

function issub(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issub(lh, label1, label2)
end

function super(s::AbstractSystem, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _super(lh, label)
end

function sub(s::AbstractSystem, hname::AbstractString, label::HLabel)
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
