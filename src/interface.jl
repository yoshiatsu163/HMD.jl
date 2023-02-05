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

function add_bond!(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer; bond_order::Rational=1//1)
    if !(0//1 <= bond_order <= 8//1)
        error("bond order must be in [0, 8] ")
    end
    topo = topology(s)
    @assert add_edge!(topo, atom_id1, atom_id2, bond_order)

    return nothing
end

function bond_order(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer)
    topo = topology(s)
    if !has_edge(topo, atom_id1, atom_id2)
        error("There is no bond beteen atoms $(atom_id1) and $(atom_id2). ")
    end

    return get_weight(topo, atom_id1, atom_id2)
end

function valence(s::AbstractSystem, atom_id::Integer)
    topo = topology(s)
    valence = 0//1
    for neigh_id in all_neighbors(topo, atom_id)
        valence += get_weight(topo, atom_id, neight_id)
    end

    return valence
end

function atom_label(atom_id::Integer)
    return HLabel("", atom_id)
end

function is_atom(label::HLabel)
    return type(label) == Category{HLabel}("")
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

function hmdsave(name::AbstractString, s::AbstractSystem)
    jldsave(name; system=s)
    return nothing
end

function hmdread(name::AbstractString)
    jldopen(name, "r+") do file
        file["system"]
    end |> return
end
