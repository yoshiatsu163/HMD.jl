# presetとしてpolymer hierarchyを追加: polymer >: {monomer >: connection, endcap >: endatom}
# add_label, add, \oplus, add!を作成 テスト可能
# selectionなど読み込み機能を追加
#
#

module DataTypes

using DataStructures
using Graphs
using LinearAlgebra
using MLStyle
using PeriodicTable
using SimpleWeightedGraphs
using StaticArrays

import Base: getindex, setproperty!, iterate, length, precision
import Base: >, <, >=, <=, +, -, *, /, ==, string, show, convert
import Base: position, time, contains
import Base: promote_rule, promote_type
import Base: ∈, ∉

include("util.jl")
include("HierarchyLabels/HierarchyLabels.jl")

using  .HierarchyLabels

# core subtype signature
export Position, BoundingBox, HLabel, H_Label, LabelHierarchy
export >, <, >=, <=, +, -, *, /, ==, string, show, convert, getindex, convert, iterate
export id, type, ==, promote_rule, promote_type, length

# core immut signature
export AbstractSystemType, GeneralSystem, AbstractSystem, System
export contains
export all_elements, element
export time, natom, nbond, topology, box, dimension, precision, system_type, show
export all_positions, position, all_travels, travel, wrapped
export hierarchy_names, hierarchy
export all_labels, count_label, has_relation, issuper, issub, super, sub, print_to_string
export prop_names, props, prop, labels_in_prop
export PackedHierarchy, rconvert, wconvert, SerializedElement

# core mut signature
export set_time!, set_box!
export _add_element!, _add_elements!, set_element!, _add_position!, _add_positions!, set_position!, set_travel!, _change_wrap!
export add_hierarchy!, remove_hierarchy!
export add_prop!, set_prop!
export add_label!, add_labels!, add_relation!, insert_relation!, add_relations!, remove_label!, remove_relation!

# trajectory specific signature


# fileIO
export SerializedTopology, PackedHierarchy, serialize, deserialize

#constants
export Entire_System, BO_Precision

const Entire_System = HLabel("entire_system", 1)
const BO_Precision = Int8

#####
##### Type `Position` definition
#####

const Position{D, F} = Vector{SVector{D, F}} where {D, F <: AbstractFloat}

function Position{D, F}(natom::Integer) where {D, F}
    return [SVector{D, F}(zeros(F, D)) for _ in 1:natom]
end

function Position{D, F}() where {D, F}
    return SVector{D, F}[]
end

function Position{D, F}(x::Matrix{F}) where {D, F}
    if size(x, 1) != D
        error("Dimension mismatch")
    end
    return [SVector{D, F}(x[:, i]) for i in 1:size(x, 2)]
end

function Position()
    Position{3, Float64}()
end

function Position(n::Integer)
    Position{3, Float64}(n)
end


#####
##### Type `BoundingBox` definition
#####

struct BoundingBox{D, F <: AbstractFloat}
    origin::SVector{D, F}
    axis::SMatrix{D, D, F}
end

function BoundingBox{D, F}(origin::Vector{F}, axis::Matrix{F}) where {D, F<:AbstractFloat}
    d = det(axis)
    if d < 0
        error("left-handed system not allowed")
    elseif d == 0
        error("box is degenerated")
    end
    return BoundingBox{D, F}(SVector{D, F}(origin), SMatrix{D, D, F}(axis))
end

function BoundingBox{D, F}() where {D, F<:AbstractFloat}
    BoundingBox{D, F}(zeros(F, 3), Matrix{F}(I, D, D))
end

#####
##### Type `System` definition
#####

abstract type AbstractSystem{D, F<:AbstractFloat} end
abstract type AbstractSystemType end

struct GeneralSystem <: AbstractSystemType end

mutable struct System{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractSystem{D, F}
    time::F
    topology::SimpleWeightedGraph{Int64, Rational{BO_Precision}}
    box::BoundingBox{D, F}

    # atom property
    position::Position{D, F}
    travel::Vector{SVector{D, Int16}}
    wrapped::Bool
    element::Vector{String}

    hierarchy::Dict{String, LabelHierarchy}
    props::Dict{String, Dict{HLabel, Any}}
end

include("property.jl")
include("trajectory.jl")

function System{D, F, SysType}() where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    System{D, F, SysType}(
        zero(F),
        SimpleWeightedGraph{Int64, Rational{BO_Precision}}(),
        BoundingBox{D, F}(),
        Position{D, F}(),
        Vector{SVector{D, Int16}}(undef, 0),
        false,
        String[],
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

function precision(s::AbstractSystem{D, F}) where {D, F<:AbstractFloat}
    return F
end

function system_type(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return SysType
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
    _add_element!(s, ename)
end

function _add_elements!(s::AbstractSystem, enames::AbstractVector{<:AbstractString})
    append!(s.element, enames)
end

function set_element!(s::AbstractSystem, atom_id::Integer, ename::AbstractString)
    s.element[atom_id] = ename
end

function set_elements!(s::AbstractSystem, atom_ids::AbstractVector{<:Integer}, enames::AbstractVector{<:AbstractString})
    if length(atom_ids) != length(enames)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.element[atom_ids] .= enames
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

    return nothing
end

function _add_positions!(s::AbstractSystem, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    append!(all_positions(s), x)
    if wrapped(s)
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    append!(s.travel, [zeros(Int16, 3) for _ in 1:length(x)])

    return nothing
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

function all_travels(s::AbstractSystem)
    return s.travel
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
    s.hierarchy |> keys |> collect
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

function add_label!(s::AbstractSystem, hname::AbstractString, label_type::AbstractString)
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

function count_label(s::AbstractSystem, hname::AbstractString, label_type::String)
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

# data conversion for HDF5

function serialize(vec::Vector{SVector{D, T}}) where {D, T<:Real}
    #return [vec[atom_id][dim] for dim in 1:D, atom_id in eachindex(vec)]
    spos = Matrix{T}(undef, D, length(vec))
    for atom_id in eachindex(vec)
        spos[:, atom_id] = vec[atom_id]
    end
    return spos
end

function deserialize(D::Integer, mat::Matrix{T}) where {T<:Real}
    pos = SVector{D, T}[]
    resize!(pos, size(mat, 2))
    for atom_id in eachindex(pos)
        pos[atom_id] = mat[:, atom_id]
    end
    return pos
end

struct SerializedTopology
    num_node::Int64
    edges_org::Vector{Int64}
    edges_dst::Vector{Int64}
    denominator::Vector{BO_Precision}
    numerator::Vector{BO_Precision}
end

function serialize(topo::SimpleWeightedGraph)
    num_node = nv(topo)
    edges_org = Vector{Int64}(undef, ne(topo))
    edges_dst = Vector{Int64}(undef, ne(topo))
    denominator = Vector{Int16}(undef, ne(topo))
    numerator = Vector{Int16}(undef, ne(topo))

    for (i, edge) in enumerate(edges(topo))
        edges_org[i], edges_dst[i] = src(edge), dst(edge)
        weight = get_weight(topo, edges_org[i], edges_dst[i])
        denominator[i], numerator[i] = weight.den, weight.num
    end

    return SerializedTopology(num_node, edges_org, edges_dst, denominator, numerator)
end

function deserialize(ser_topo::SerializedTopology)
    topo = SimpleWeightedGraph{Int64, Rational{BO_Precision}}()
    add_vertices!(topo, ser_topo.num_node)
    for i in 1:length(ser_topo.edges_org)
        add_edge!(topo, ser_topo.edges_org[i], ser_topo.edges_dst[i], Rational{BO_Precision}(ser_topo.numerator[i], ser_topo.denominator[i]))
        #add_edge!(topo, ser_topo.edges_org[i], ser_topo.edges_dst[i], ser_topo.numerator[i]//ser_topo.denominator[i])
    end

    return topo
end

include("test.jl")

"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""



end #module
