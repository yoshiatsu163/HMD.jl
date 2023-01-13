# heierchy最下位にatomをいれると辞書をなくせる
# Dict一発呼び出しよりも速度は遅くなる
# 以下コードはキャッシュ用途に保管？
mutable struct Label2atom{T1, T2}
    toatom::Dict{T1, Set{T2}}
    tolabel::Dict{T2, Set{T1}}
end

function Label2atom()
    Label2atom{Label, Int64}(
        Dict{Label, Set{Int64}}(),
        Dict{Int64, Set{Label}}()
    )
end

function toatom(l2a::Label2atom)
    l2a.toatom
end

function tolabel(l2a::Label2atom)
    l2a.tolabel
end

function _labels(l2a::Label2atom)
    toatom(l2a) |> keys |> Set
end

function _atoms(l2a::Label2atom)
    tolabel(l2a) |> keys |> Set
end

function getindex(l2a::Label2atom, label::Label)
    toatom(l2a)[label]
end

function getindex(l2a::Label2atom, atomid::Integer)
    tolabel(l2a)[atomid]
end

# add labelのみだとlabelのついてない原子を扱えない
# 1. Heierchyにそもそも入れない
# 2. add_atom等を用意してatomとのコヒーレンシを保つ
#   sync!()でもいいがテストユニットが大きくなる
function _add_label!(l2a::Label2atom, label::Label, ids::Vector{<:Integer})
    if label ∈ _labels(l2a)
        error("label already exists in Label => atom mapping. ")
    end
    nf = [id for id in ids if id ∉ _atoms(l2a)]
    if !isempty(nf)
        str = join(string.(nf), " ")
        error("atom $str not found in Label => atom mapping. ")
    end

    push!(toatom(l2a) , label => Set(ids))
    for id in ids
        push!(tolabel(l2a)[id], label)
    end
end

function _add_atoms!(l2a::Label2atom, n::Integer)
    Natom = tolabel(l2a) |> length
    for i in 1:n
        push!(tolabel(l2a), Natom+i => Set{Label}())
    end
    tolabel(l2a) |> length
end

function _add_atom!(l2a::Label2atom)
    _add_atoms!(l2a, 1)
end

function ==(lhs::Label2atom, rhs::Label2atom)
    toatom(lhs) == toatom(rhs) && tolabel(lhs) == tolabel(rhs)
end
