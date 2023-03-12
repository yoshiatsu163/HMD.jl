# openmm trajectoryもバックエンドにしたいのでDataTypesまわりの最小apiを決める必要がある
const Entire_System = HLabel("entire_system", 1)

const atom_mass = Dict{StaticString, Float64}(
    StaticString(elements[:H ].symbol) => 1.008,
    StaticString(elements[:C ].symbol) => 12.012,
    StaticString(elements[:N ].symbol) => 14.007,
    StaticString(elements[:O ].symbol) => 16.000,
    StaticString(elements[:F ].symbol) => 19.000,
    StaticString(elements[:Si].symbol) => 28.086,
    StaticString(elements[:S ].symbol) => 32.067,
    StaticString(elements[:Cl].symbol) => 35.453
)

function dimension(s::AbstractSystem{D, F}) where {D, F<:AbstractFloat}
    return D
end

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::AbstractString; super::HLabel)
    add_atom!(s, x, elem; super=super)

    return nothing
end

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::Category{Element}; super::HLabel)
    atom_id = natom(s) + 1
    _add_position!(s, x)
    _add_element!(s, elem)
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        add_label!(s, hname, atom_label(atom_id))
        add_relation!(s, hname; super=super, sub=atom_label(atom_id))
    end

    return nothing
end

function add_atoms!(s::AbstractSystem, x::AbstractVector{<:AbstractVector{<:AbstractFloat}}, elem::AbstractVector{<:AbstractString}; super::HLabel)
    add_atoms!(s, x, Category{Element}.(elem); super=super)

    return nothing
end

function add_atoms!(s::AbstractSystem, x::AbstractVector{<:AbstractVector{<:AbstractFloat}}, elem::AbstractVector{Category{Element}}; super::HLabel)
    if length(elem) != length(x)
        error("length of elem and x must be same. ")
    end

    front = natom(s) + 1
    back = natom(s) + length(x)
    atom_labels = atom_label.(front:back)
    for i in 1:length(x)
        _add_position!(s, x[i])
        _add_element!(s, elem[i])
    end
    @assert add_vertices!(topology(s), length(x))
    for hname in hierarchy_names(s)
        add_labels!(s, hname, atom_label.(front:back))
        add_relations!(s, hname; super=super, subs=atom_labels)
    end

    return nothing
end

function add_bond!(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer; bond_order::Rational=1//1)
    if !(0//1 <= bond_order <= 8//1)
        error("bond order must be in [0, 8] ")
    end
    topo = topology(s)
    @assert add_edge!(topo, atom_id1, atom_id2, bond_order)

    return nothing
end

#function add_bonds!(s::AbstractSystem, pair::AbstractVector{Tuple{<:Integer, <:Integer, Rational{<:Integer}}})
function add_bonds!(s::AbstractSystem, pair)
    if any(p -> !(0//1 <= p[3] <= 8//1), pair)
        error("bond order must be in [0, 8] ")
    end

    topo = topology(s)
    for (atom_id1, atom_id2, bond_order) in pair
        @assert add_edge!(topo, atom_id1, atom_id2, bond_order)
    end
end

function add_bond!(s::AbstractSystem, label1::HLabel, label2::HLabel; bond_order::Rational=1//1)
    if !is_atom(label1) || !is_atom(label2)
        error("label is not for atom. ")
    end

    atom_id1 = convert(Int64, id(label1))
    atom_id2 = convert(Int64, id(label2))
    add_bond!(s, atom_id1, atom_id2; bond_order=bond_order)

    return nothing
end

function bond_order(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer)
    topo = topology(s)
    if !has_edge(topo, atom_id1, atom_id2)
        error("There is no bond beteen atoms $(atom_id1) and $(atom_id2). ")
    end

    return get_weight(topo, atom_id1, atom_id2)
end

function bond_order(s::AbstractSystem, label1::HLabel, label2::HLabel)
    if !is_atom(label1) || !is_atom(label2)
        error("label is not for atom. ")
    end
    atom_id1 = convert(Int64, id(label1))
    atom_id2 = convert(Int64, id(label2))

    return bond_order(s, atom_id1, atom_id2)
end

function valence(s::AbstractSystem, atom_id::Integer)
    topo = topology(s)
    valence = 0//1
    for neigh_id in all_neighbors(topo, atom_id)
        valence += get_weight(topo, atom_id, neigh_id)
    end

    return valence
end

function valence(s::AbstractSystem, label::HLabel)
    if !is_atom(label)
        error("label $label is not for atom. ")
    end

    return valence(s, convert(Int64, id(label1)))
end

function atom_label(atom_id::Integer)
    return HLabel("", atom_id)
end

function is_atom(label::HLabel)
    return type(label) == Category{H_Label}("")
end

function neighbors(s::AbstractSystem, atom_id::Integer)
    topo = topology(s)
    return all_neighbors(topo, atom_id)
end

function neighbors(s::AbstractSystem, label::HLabel)
    if !is_atom(label)
        error("label $label is not for atom. ")
    end
    topo = topology(s)

    return all_neighbors(topo, convert(Int64, id(label1)))
end

function l2a(s::AbstractSystem, hname::AbstractString, label::HLabel)
    sub_labels = [label]
    atom_ids = Vector{Int64}(undef, 0)

    while any(!is_atom, sub_labels)
        sub_next = Vector{HLabel}(undef, 0)
        for l in sub_labels
            if is_atom(l)
                push!(atom_ids, convert(Int64, id(l)))
            else
                append!(sub_next, sub(s, hname, l))
            end
        end
        sub_labels = sub_next
    end
    append!(atom_ids, convert.(Int64, id.(sub_labels)))

    return atom_ids
end

# TODO
#function atom_label(T::Type{<:AbstractLabel}, atom_id::Integer)
#    T(atom_id, "")
#end

#function add!(s::System, addend::System)
    #topology の connection pointが必要
    #一般のsystemは自由度が高すぎてaddが定義困難 -> HMDPolymer等を作って限定的にaddを定義
    #   AbstractSystemの定義が必要
#end

# 同名の関数がDataTypesにあり
function all_labels(s::AbstractSystem, hname::AbstractString, label_type::Category{H_Label})
    labels = hierarchy(s, hname) |> DataTypes.HierarchyLabels._label2node |> keys |> collect

    return filter!(label -> type(label)==label_type, labels)
end

function super_labels(s::AbstractSystem, hname::AbstractString, label::HLabel)
    return _traverse_from(s, hname, label, super)
end

function sub_labels(s::AbstractSystem, hname::AbstractString, label::HLabel)
    return _traverse_from(s, hname, label, sub)
end

function _traverse_from(s::AbstractSystem, hname::AbstractString, label::HLabel, func::Function)
    labels = Vector{HLabel}(undef, 0)

    # func is either super or sub
    stack = [label]
    while !isempty(stack)
        current = popfirst!(stack)
        next = func(s, hname, current)
        prepend!(stack, next)
        push!(labels, current)
    end

    return labels
end

function wrap!(s::AbstractSystem)
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

function unwrap!(s::AbstractSystem)
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

function hmdsave(name::AbstractString, s::AbstractSystem{D, F}; compress=false) where {D, F<:AbstractFloat}
    if compress
        println("warning: compression is not supported yet.")
    end

    h5open(name, "w") do file
        # metadata as type D, F, Systype
        file["time"] = time(s)
        file["position"] = serialize(all_positions(s))
        file["travel"] = serialize(all_travels(s))
        file["wrapped"] = wrapped(s)

        file["box/origin"] = Vector(box(s).origin)
        file["box/axis"] = Matrix(box(s).axis)

        selem = serialize(all_elements(s))
        file["element/chars"]  = selem.chars
        file["element/bounds"] = selem.bounds

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
    end

    return nothing
end

function hmdread(name::AbstractString)
    #return jldopen(name, "r") do file
    #    system_type = file["system_type"]
    #    s = system_type()
    #    D, F = dimension(s), precision(s)
#
    #    set_time!(s, file["time"])
    #    set_box!(s, file["box"])
    #    s.position = rconvert(Position{D, F}, file["position"])
    #    s.travel = rconvert(Vector{MVector{D, Int16}}, file["travel"])
    #    s.wrapped = file["wrapped"]
    #    s.element = rconvert(Vector{Category{Element}}, file["element"])
    #    s.topology = file["topology"]
    #    for hname in file["hierarchy_names"]
    #        add_hierarchy!(s, hname)
    #        s.hierarchy[hname] = rconvert(LabelHierarchy, file["hierarchy/$hname"])
    #    end
    #    s.props = file["props"]
    #    s
    #end
end
