const Entire_System = Label(1, "entire_system")

function add_atom!(s::System, x::AbstractVector{<:AbstractFloat}, elem::AbstractString)
    atom_id = natom(s) + 1
    push!(s.position, x)
    push!(s.element , Category{Element}(elem))
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        add_label!(s, hname, atom_label(atom_id))
        add_relation!(s, hname; super=Entire_System, sub=atom_label(atom_id))
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

function add!(s::System, addend::System)
    #topology の connection pointが必要
    #一般のsystemは自由度が高すぎてaddが定義困難 -> HMDPolymer等を作って限定的にaddを定義
    #   AbstractSystemの定義が必要
end

#function ⊕(lhs::System, rhs::System)
#    # 重複あればエラー
#    # 同一名のHierarchyあればadd
#    # 異なるHierarchyあれば追加
#end

#function ⊗=(lhs::System, rhs::System)
#    add!(s::System, addend::System)
#    return nothing
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
