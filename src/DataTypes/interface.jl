# IOメモ
# 行のフォーマットを指定してdatatypesに読み書きできるようにする
# IOメモ
# IOメモ

function add_atom!(s::System, x::AbstractVector{<:AbstractFloat}, elem::AbstractString)
    atom_id = natom(s) + 1
    push!(s.position, x)
    push!(s.element , Category{Element}(elem))
    @assert add_vertex!(topology(s))
    for hname in hierarchy_names(s)
        lh = hierarchy(s, hname)
        _add_label!(lh, atom_label(atom_id))
        _add_relation!(lh; super=Entire_System, sub=atom_label(atom_id))
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
