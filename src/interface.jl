const Entire_System = HLabel("entire_system", 1)

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::AbstractString; super::HLabel)
    atom_id = natom(s) + 1
    push!(s.position, x)
    push!(s.element , Category{Element}(elem))
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        add_label!(s, hname, atom_label(atom_id))
        add_relation!(s, hname; super=super, sub=atom_label(atom_id))
    end

    return nothing
end

function add_atom!(s::AbstractSystem, x::AbstractVector{<:AbstractFloat}, elem::Category{Element}; super::HLabel)
    add_atom!(s, x, string(elem); super=super)
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
    return type(label) == Category{HLabel}("")
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
function all_labels(s::AbstractSystem, hname::AbstractString, label_type::Category{HLabel})
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

@inline function wrap(box::BoundingBox, x::AbstractVector{<:AbstractFloat})
    wrapped = x .- box.origin

    coeff = [cdot(wrapped, box.axis[:, i]) for i in 1:length(x)] # x = reduce(+, coeff[i] .* axis[:, i] for i in 1:length(x))
    coeff .-= floor.(coeff)

    for dim in 1:length(x)
        wrapped .+= coeff[dim] .* box.axis[:, dim]
    end

    return box.origin .+ wrapped
end

function hmdsave(name::AbstractString, s::AbstractSystem)
    jldsave(name; system=s)
    return nothing
end

function hmdread(name::AbstractString)
    jldopen(name, "r+") do file
        file["system"]
    end |> return
end
