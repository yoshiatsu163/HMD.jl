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
using SimpleWeightedGraphs
using StaticArrays

import Base: getindex, setproperty!, iterate, length, precision, similar
import Base: >, <, >=, <=, +, -, *, /, ==, string, show, convert
import Base: position, time, contains
import Base: promote_rule, promote_type
import Base: ∈, ∉
import Base: close

include("util.jl")
include("HierarchyLabels/HierarchyLabels.jl")

using  .HierarchyLabels

# core subtype signature
export Position, BoundingBox, HLabel, H_Label, LabelHierarchy
export >, <, >=, <=, +, -, *, /, ==, string, show, convert, getindex, convert, iterate
export id, type, ==, promote_rule, promote_type, length, similar

# core immut signature
export AbstractSystemType, GeneralSystem, AbstractSystem, System
export contains
export all_elements, element
export time, natom, nbond, topology, box, dimension, precision, system_type, show
export all_positions, position, all_travels, travel, wrapped
export hierarchy_names, hierarchy
export all_labels, count_label, has_relation, issuper, issub, super, sub, print_to_string
export prop_names, props, prop, labels_in_prop

# core mut signature
export set_time!, set_box!
export _add_position!, _add_positions!, set_position!, set_travel!, _change_wrap!
export _add_element!, _add_elements!, set_element!, set_elements!
export add_hierarchy!, remove_hierarchy!
export add_prop!, set_prop!
export add_label!, add_labels!, add_relation!, insert_relation!, add_relations!, remove_label!, remove_relation!

# trajectory specific signature


# fileIO
export AbstractH5, H5system, H5traj, SerializedTopology, PackedHierarchy
export serialize, deserialize, h5system, h5traj, get_file, close

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
    @assert natm == length(s.element)
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

function hierarchy_names(s::System)
    s.hierarchy |> keys |> collect
end

function hierarchy(s::System, hname::AbstractString)
    s.hierarchy[hname]
end

function add_hierarchy!(s::System, hname::AbstractString)
    if hname in hierarchy_names(s)
        error("hierarchy $(hname) already exists. ")
    end
    push!(s.hierarchy, hname => LabelHierarchy())
    #_add_label!(hierarchy(s, hname), Entire_System)
    add_label!(s, hname, Entire_System)

    return nothing
end

function remove_hierarchy!(s::System, hname::AbstractString)
    delete!(s.hierarchy, hname)
end

function all_labels(s::System, hname::AbstractString)
    lh = hierarchy(s, hname)
    return _labels(lh)
end

function all_labels(s::AbstractSystem, hname::AbstractString, label_type::AbstractString)
    labels = hierarchy(s, hname) |> DataTypes.HierarchyLabels._label2node |> keys |> collect

    return filter!(label -> type(label)==label_type, labels)
end

function add_label!(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    result = _add_label!(lh, label)
    @match result begin
        Label_Occupied => error("label $(label) already exists. ")
        Success        => return nothing
        _              => error("fatal error")
    end
end

function add_label!(s::System, hname::AbstractString, label_type::AbstractString)
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

function add_labels!(s::System, hname::AbstractString, label_types::AbstractVector{<:HLabel})
    lh = hierarchy(s, hname)

    _add_labels!(lh, label_types)

    return nothing
end

function count_label(s::System, hname::AbstractString, label_type::String)
    lh = hierarchy(s, hname)
    labels = _label2node(lh) |> keys
    return count(l -> type(l)==label_type, labels)
end

function add_relation!(s::System, hname::AbstractString; super::HLabel, sub::HLabel)
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

function add_relations!(s::System, hname::AbstractString; super::HLabel, subs::AbstractVector{HLabel})
    lh = hierarchy(s, hname)
    for sub in subs
        _add_relation!(lh; super=super, sub=sub, unsafe=true)
    end

    return nothing
end

function insert_relation!(s::System, hname::AbstractString, label::HLabel; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    add_label!(s, hname, label)
    add_relation!(s, hname; super=super, sub=label)
    add_relation!(s, hname; super=label, sub=sub)
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function remove_label!(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)

    if !_contains!(lh, label)
        error("label $(label) not found. ")
    end
    @assert _remove_label!(lh, label)

    return nothing
end

function remove_relation!(s::System, hname::AbstractString; super::HLabel, sub::HLabel)
    lh = hierarchy(s, hname)

    if !_has_relation(lh, super, sub)
        error("relation betewwn $(super) and $(sub) not found. ")
    end
    @assert _remove_relation!(lh, super, sub)

    return nothing
end

function contains(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _contains(lh, label)
end

function issuper(s::System, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issuper(lh, label1, label2)
end

function issub(s::System, hname::AbstractString, label1::HLabel, label2::HLabel)
    lh = hierarchy(s, hname)
    return _issub(lh, label1, label2)
end

function super(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _super(lh, label)
end

function sub(s::System, hname::AbstractString, label::HLabel)
    lh = hierarchy(s, hname)
    return _sub(lh, label)
end

function Base.similar(s::System{D, F, SysType}) where {D, F, SysType}
    return System{D, F, SysType}()
end

#####
##### HDF5 IO
#####

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

abstract type AbstractH5 end

mutable struct H5system <: AbstractH5
    file::Union{HDF5.File, HDF5.Group}
end

mutable struct H5traj <: AbstractH5
    file::Union{HDF5.File, HDF5.Group}
end

function h5system(name::AbstractString, mode::AbstractString)
    file_handler = H5system(h5open(name, mode))
    file = get_file(file_handler)
    if mode != "w" && read(file, "infotype") != "System"
        close(file_handler)
        error("file $(name) is not a System file. ")
    end

    return file_handler
end

function h5traj(name::AbstractString, mode::AbstractString)
    file_handler = H5traj(h5open(name, mode))
    file = get_file(file_handler)
    if mode != "w" && read(file, "infotype") != "Trajectory"
        close(file_handler)
        error("file $(name) is not a Trajectory file. ")
    end

    return file_handler
end

function close(file_handler::AbstractH5)
    close(get_file(file_handler))
end

function get_file(file_handler::AbstractH5)
    return file_handler.file
end

function hmdsave(file_handler::H5system, s::System{D, F, SysType}; compress=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if compress
        println("warning: compression is not supported yet.")
    end

    file = get_file(file_handler)

    # metadata
    file["infotype"] = "System"
    file["dimension"] = dimension(s)
    file["precision"] = precision(s) |> string
    file["system_type"] = system_type(s) |> string

    #data
    file["time"] = time(s)
    file["position"] = serialize(all_positions(s))
    file["travel"] = serialize(all_travels(s))
    file["wrapped"] = wrapped(s)

    file["box/origin"] = Vector(box(s).origin)
    file["box/axis"] = Matrix(box(s).axis)

    chars, bounds = serialize(all_elements(s))
    file["element/chars"]  = chars
    file["element/bounds"] = bounds

    stopo = serialize(topology(s))
    for fname in fieldnames(SerializedTopology)
        file["topology/$(fname)"] = getfield(stopo, fname)
    end

    file["hierarchy_names"] = hierarchy_names(s)
    for hname in hierarchy_names(s)
        ser_hierarchy = serialize(hierarchy(s, hname))
        for fname in fieldnames(PackedHierarchy)
            file["hierarchy/$hname/$(fname)"] = getfield(ser_hierarchy, fname)
        end
    end
    ##temporary
    #file["props"] = s.props

    return nothing
end

function read_system(system_file::H5system, template::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    # metadata check
    D_file, F_file, SysType_file = get_metadata(system_file)
    if (D_file, F_file) != (D, string(F))
        error("System type mismatch. (file: $(D_file), $(F_file), $(SysType_file), system: $(D), $(F), $(SysType)")
    end

    # data
    s = similar(template)
    import_dynamic!(s, system_file)
    import_static!(s, system_file)
    return s
end

function import_dynamic!(s::System{D, F, SysType}, system_file::H5system) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    file = get_file(system_file)
    set_time!(s, read(file, "time"))
    set_box!(s, BoundingBox{D, F}(read(file, "box/origin"), read(file, "box/axis")))
    s.position = deserialize(D, read(file, "position"))
    s.travel = deserialize(D, read(file, "travel"))

    return nothing
end

function import_static!(s::System{D, F, SysType}, system_file::H5system) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    file = get_file(system_file)
    s.element = deserialize(read(file, "element/chars"), read(file, "element/bounds"))
    s.wrapped = read(file, "wrapped")
    s.topology = SerializedTopology(read(file, "topology/num_node"),
                                    read(file, "topology/edges_org"),
                                    read(file, "topology/edges_dst"),
                                    read(file, "topology/denominator"),
                                    read(file, "topology/numerator")) |> deserialize
    for hname in read(file, "hierarchy_names")
        if hname ∉ hierarchy_names(s)
            add_hierarchy!(s, hname)
        end
        ph = PackedHierarchy(read(file, "hierarchy/$hname/num_node"),
                            read(file, "hierarchy/$hname/edges_org"),
                            read(file, "hierarchy/$hname/edges_dst"),
                            read(file, "hierarchy/$hname/label_ids"),
                            read(file, "hierarchy/$hname/chars"),
                            read(file, "hierarchy/$hname/bounds"))
        s.hierarchy[hname] = deserialize(ph)
    end

    return nothing
end

function get_metadata(system_file::H5system)
    file = get_file(system_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type")

    return D, F, SysType
end

include("trajectory.jl")
include("property.jl")
include("test.jl")

"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""



end #module
