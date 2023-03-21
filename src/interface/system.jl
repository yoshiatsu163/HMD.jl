# openmm trajectoryもバックエンドにしたいのでDataTypesまわりの最小apiを決める必要がある
const Entire_System = HLabel("entire_system", 1)

const atom_mass = Dict{String, Float64}(
    elements[:H ].symbol => 1.008,
    elements[:C ].symbol => 12.012,
    elements[:N ].symbol => 14.007,
    elements[:O ].symbol => 16.000,
    elements[:F ].symbol => 19.000,
    elements[:Si].symbol => 28.086,
    elements[:S ].symbol => 32.067,
    elements[:Cl].symbol => 35.453
)

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::AbstractString; super::HLabel)
    atom_id = natom(s) + 1
    DataTypes._add_position!(s, x)
    DataTypes._add_element!(s, elem)
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        add_label!(s, hname, atom_label(atom_id))
        add_relation!(s, hname; super=super, sub=atom_label(atom_id))
    end

    return nothing
end

function add_atoms!(s::AbstractSystem, x::AbstractVector{<:AbstractVector{<:AbstractFloat}}, elem::AbstractVector{<:AbstractString}; super::HLabel)
    if length(elem) != length(x)
        error("length of elem and x must be same. ")
    end

    front = natom(s) + 1
    back = natom(s) + length(x)
    atom_labels = atom_label.(front:back)
    DataTypes._add_positions!(s, x)
    DataTypes._add_elements!(s, elem)
    @assert add_vertices!(topology(s), length(x))
    for hname in hierarchy_names(s)
        add_labels!(s, hname, atom_labels)
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
    return type(label) == ""
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

#####
##### System HDF5 interface
#####

function hmdsave(name::AbstractString, s::AbstractSystem{D, F}; compress=false) where {D, F<:AbstractFloat}
    file = h5system(name, "w")
    DataTypes.hmdsave(file, s; compress=compress)
    close(file)

    return nothing
end

function hmdread!(s::AbstractSystem{D, F}, name::AbstractString) where {D, F<:AbstractFloat}
    file = h5system(name, "r")
    DataTypes.read_system(file)
    close(file)

    return s
end
 
end