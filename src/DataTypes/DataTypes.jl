# presetとしてpolymer hierarchyを追加: polymer >: {monomer >: connection, endcap >: endatom}
# add_label, add, \oplus, add!を作成 テスト可能
# selectionなど読み込み機能を追加
#
#

# ここにはデータ型の内部を直接触る最小限の関数を書く
# HMD/interface/system.jlにAbstractSystem関数を作成することでDataTypesの内部構造を容易に変更できる

module DataTypes

using DataStructures
using Graphs
using HDF5
using LinearAlgebra
using MLStyle
using PeriodicTable
using Reexport
using SimpleWeightedGraphs
using StaticArrays

using ..HierarchyLabels

@reexport import Base: *, +, -, /, <, <=, ==, >, >=, close, contains, convert, getindex,firstindex, lastindex, iterate,
    length, position, precision, promote_rule, promote_type, setproperty!, show, similar,
    string, time, ∈, ∉

@reexport import ..HMD: deserialize, serialize
@reexport import ..HMD:
    # system core interface
    AbstractBbox,
    AbstractSystem,
    AbstractSystemType,
    AbstractTrajectory,
    dimension,
    precision,
    system_type,
    similar,
    show,
    natom,
    nbond,
    time,
    set_time!,
    topology,
    box,
    set_box!,
    all_elements,
    element,
    element,
    set_element!,
    set_elements!,
    all_positions,
    position,
    set_position!,
    set_position!,
    all_travels,
    travel,
    set_travel!,
    wrapped,
    wrap!,
    unwrap!,
    label2atom,

    # system label manipulation
    hierarchy_names,
    hierarchy,
    add_hierarchy!,
    remove_hierarchy!,
    all_labels,
    add_label!,
    add_labels!,
    count_label,
    add_relation!,
    add_relations!,
    insert_relations!,
    remove_label!,
    remove_relation!,
    contains,
    issuper,
    issub,
    super,
    sub,

    # system io interface
    AbstractFileFormat,
    close,

    # trajectory interface
    get_system,
    all_timesteps,
    get_timestep,
    length,
    add!,
    import_dynamic!,
    import_static!,
    latest_reaction,
    similar,
    similar_system,
    dimension,
    precision,
    system_type,
    wrapped,
    wrap!,
    unwrap!,
    add_snapshot!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    getindex,
    lastindex,
    firstindex,
    wrapped,

    # trajectory io interface
    add_snapshot!,
    import_dynamic!,
    import_static!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    getindex,
    lastindex,
    firstindex,
    wrapped

# core subtype signature
export Position, BoundingBox, HLabel, LabelHierarchy

# core immut signature
export GeneralSystem, System, print_to_string

# fileIO
export H5system, H5traj, SerializedTopology, PackedHierarchy
export h5system, h5traj, get_file

#constants
export Entire_System, BO_Precision

const Entire_System = HLabel("entire_system", 1)
const BO_Precision = Int8

include("position.jl")
include("boundingbox.jl")

#####
##### Type `System` definition
#####

struct GeneralSystem <: AbstractSystemType end

mutable struct System{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractSystem{D, F, SysType}
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

function dimension(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return D
end

function precision(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return F
end

function system_type(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return SysType
end

function Base.similar(s::System{D, F, SysType}) where {D, F, SysType}
    return System{D, F, SysType}()
end

function Base.show(io::IO, ::MIME"text/plain", s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    "System{$D, $F}
        time: $(time(s))
        bbox: $(box(s))
        natoms: $(natom(s))
        hierarchy: $(hierarchy_names(s))
    " |> println
end

function natom(s::System)
    natm = length(s.position)
    #@assert natm == length(s.element)
    return length(s.position)
end

function nbond(s::System)
    return topology(s) |> ne
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

function all_elements(s::System)
    s.element
end

function element(s::System, atom_id::Integer)
    s.element[atom_id]
end

function element(s::System, label::HLabel)
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    s.element[atom_id]

    return nothing
end

function _add_element!(s::System, ename::AbstractString)
    _add_element!(s, ename)
end

function _add_elements!(s::System, enames::AbstractVector{<:AbstractString})
    append!(s.element, enames)
end

function set_element!(s::System, atom_id::Integer, ename::AbstractString)
    s.element[atom_id] = ename
end

function set_elements!(s::System, atom_ids::AbstractVector{<:Integer}, enames::AbstractVector{<:AbstractString})
    if length(atom_ids) != length(enames)
        throw(DimensionMismatch("Dimension of atom_ids is $(length(atom_ids)) but enames dimension is $(length(enames))"))
    end
    s.element[atom_ids] .= enames
end

function all_positions(s::System)
    s.position
end

function position(s::System, atom_id::Integer)
    s.position[atom_id]
end

function _add_position!(s::System, x::AbstractVector{<:AbstractFloat})
    push!(s.position, x)
    if wrapped(s)
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    push!(s.travel, zeros(Int16, 3))

    return nothing
end

function _add_positions!(s::System, x::AbstractVector{<:AbstractVector{<:AbstractFloat}})
    append!(all_positions(s), x)
    if wrapped(s)
        error("atom addition with wrapped coordinates is not supprted. ")
    end
    append!(s.travel, [zeros(Int16, 3) for _ in 1:length(x)])

    return nothing
end

function set_position!(s::System, atom_id::Integer, x::AbstractVector{<:AbstractFloat})
    s.position[atom_id] = x
end

function set_position!(s::System, label::HLabel, x::AbstractVector{<:AbstractFloat})
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    set_position!(s, label, x)

    return nothing
end

function all_travels(s::System)
    return s.travel
end

function travel(s::System, atom_id::Integer)
    return s.travel[atom_id]
end

function set_travel!(s::System, atom_id::Integer, n::AbstractVector{<:Integer})
    s.travel[atom_id] = n
end

function wrapped(s::System)
    s.wrapped
end

function _change_wrap!(s::System)
    s.wrapped = !(s.wrapped)
end

function wrap!(s::System)
    if wrapped(s)
        return nothing
    end

    axis = box(s).axis
    origin = box(s).origin
    # x = c[1] .* axis[:,1] .+ c[2] .* axis[:,2] .+ ...
    e_i_e_j = [dot(axis[:,i], axis[:,j]) for i in 1:dimension(s), j in 1:dimension(s)] |> Symmetric
    for id in 1:natom(s)
        x = position(s, id) .- origin
        c = e_i_e_j \ [dot(x, axis[:,dim]) for dim in 1:dimension(s)]
        travel = floor.(Int16, c) #trunc.(Int16, c)
        set_travel!(s, id, travel)
        digit = c .- travel
        pos = map(1:dimension(s)) do dim
            #if digit[dim] >= 0
            #    digit[dim] .* axis[:,dim]
            #else
            #    (digit[dim] + 1) .* axis[:,dim]
            #end
            digit[dim] .* axis[:,dim]
        end |> p-> reduce(.+, p)
        set_position!(s, id, pos .+ origin)
    end
    _change_wrap!(s)

    return nothing
end

function unwrap!(s::System)
    if !wrapped(s)
        return nothing
    end

    axis   = box(s).axis
    origin = box(s).origin
    for i in 1:natom(s)
        x = position(s, i) .- origin
        n = travel(s, i)
        pos = x .+ mapreduce(dim -> n[dim] .* axis[:, dim], .+, 1:dimension(s))
        set_position!(s, i, pos .+ origin)
        set_travel!(s, i, zeros(Int16, 3))
    end
    _change_wrap!(s)

    return nothing
end

function label2atom(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    labels = _labels(lh)
    atom_ids = Int64[]

    stack = [_get_nodeid(lh, label)]
    if is_atom(labels[stack[1]])
        return [stack[1]]
    end

    while !isempty(stack)
        current_node = popfirst!(stack)
        next_nodes = _sub_id(lh, current_node)
        for node in next_nodes
            label = labels[node]
            if is_atom(label)
                push!(atom_ids, id(label))
            else
                pushfirst!(stack, node)
            end
        end
    end

    return atom_ids
end

function is_atom(label::HLabel)
    return type(label) == ""
end

include("label_manipulation.jl")
include("property.jl")
include("system_io.jl")
include("trajectory.jl")
include("test.jl")

end #module
