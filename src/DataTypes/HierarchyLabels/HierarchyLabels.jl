module HierarchyLabels

using Graphs
using MetaGraphs
using MLStyle

using ..DataTypes: Id, Category

import Base: getindex
import Base: ==
import Base: ∈
import Base: ∋
import Base: ∉
import Base: show

export HLabel, LabelHierarchy
export LabelResult, Label_Missing, Label_Occupied, Label_Duplication, Relation_Missing, Relation_Occupied, Success
export id, type, ==
export _add_label!, _add_relation!, _remove_label!, _remove_relation!
export _label2node, _contains, ∈, ∋, ∉, _has_relation ,_get_nodeid, getindex
export _issuper, _issub, _super_id, _sub_id, _super, _sub
export _root, _depth

# AbstractLabel must define constructor s.t. ALabel(id::Integer, type)
abstract type AbstractLabel end

struct HLabel <: AbstractLabel
    type::Category{HLabel}
    id::Id{HLabel}
end

function HLabel(type::AbstractString, id::Integer)
    HLabel(Category{HLabel}(type), Id{HLabel}(id))
end

#const No_Label = HLabel(typemin(Int64), "No_Label")

#@enum LabelResult Label_Missing=1 Label_Occupied=2 Label_Duplication=3 Relation_Missing=4 Relation_Occupied=5 Success=6
@data LabelResult begin
    Label_Missing
    Label_Occupied
    Label_Duplication
    Relation_Missing
    Relation_Occupied
    Success
end

function id(label::HLabel)
    label.id
end

function type(label::HLabel)
    label.type
end

function ==(lhs::HLabel, rhs::HLabel)
    id(lhs) == id(rhs) && type(lhs) == type(rhs)
end

function Base.show(io::IO, ::MIME"text/plain", label::HLabel)
    i = convert(Int64, id(label))
    t = type(label) |> string
    print(io, "HLabel(\"$t\", $i)")
end

function Base.show(label::HLabel)
    i = convert(Int64, id(label))
    t = type(label) |> string
    print("HLabel(\"$t\", $i)")
end

function Base.print_to_string(label::HLabel)
    i = convert(Int64, id(label))
    t = type(label) |> string
    return "HLabel(\"$t\", $i)"
end

Base.@kwdef mutable struct LabelHierarchy
    mg::MetaDiGraph = MetaDiGraph()
    label2node::Dict{HLabel, Int64} = Dict{HLabel, Int64}()
end

function Base.show(io::IO, ::MIME"text/plain", lh::LabelHierarchy)
    subtree(lh, label, level, indent) = begin
        label == _root(lh) && println("$label")
        for s in _sub(lh, label)
            str = join(fill(" ", level * indent)) * "$s"
            println(io, str)
            subtree(lh, s, level + 1, indent)
        end
    end
    subtree(lh, _root(lh), 1, 4)
end

function _hierarchy(lh::LabelHierarchy)
    return lh.mg
end

function _label2node(lh::LabelHierarchy)
    return lh.label2node
end

function update_l2n!(lh::LabelHierarchy)
    mg = _hierarchy(lh)
    l2n = _label2node(lh)
    labels = [props(mg, i)[:label] for i in vertices(mg)]

    @assert allunique(labels)
    @assert !is_cyclic(mg)

    for (node, label) in zip(vertices(mg), labels)
        @assert label ∈ keys(l2n)
        l2n[label] = node
    end

    return nothing
end

function _contains(lh::LabelHierarchy, label::HLabel)
    #println(label)
    #println("\n")
    #for p in _label2node(lh) |> pairs
    #    println(p)
    #end
    #return label ∈ (_label2node(lh) |> keys |> collect)

    mg = _hierarchy(lh)
    for i in vertices(mg)
        label == get_prop(mg, i, :label) && return true
    end

    return false
end

function ∈(label::HLabel, lh::LabelHierarchy)
    return _contains(lh, label)
end

function ∋(lh::LabelHierarchy, label::HLabel)
    return _contains(lh, label)
end

function ∉(label::HLabel, lh::LabelHierarchy)
    return !_contains(lh, label)
end

function _get_nodeid(lh::LabelHierarchy, label::HLabel)
    if length(_hierarchy(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return _label2node(lh)[label]
end

function _get_label(lh::LabelHierarchy, id::Integer)
    if length(_hierarchy(lh)) == 0
        error("There is no label in LabelHierarchy. ")
    end
    return get_prop(_hierarchy(lh), id, :label)
end

function _labels(lh::LabelHierarchy)
    @assert nv(_hierarchy(lh)) == _label2node(lh) |> length
    return [_get_label(lh, i) for i in 1:nv(_hierarchy(lh))]
end

# TODO: ここだけ見てすべての場合について挙動がわかるようにする
function _has_relation(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !(_contains(lh, label1) && _contains(lh, label2))
        error("labels not found in LabelHierarchy. ")Label_Occupied
    end
    return _issuper(lh, label1, label2) || _issub(lh, label1, label2)
end

# ここでinsertionを定義するのは複雑なのでdatatypesのほうがよさそう
function _add_label!(lh::LabelHierarchy, label::HLabel)
    mg = _hierarchy(lh)

    if _contains(lh, label)
        return Label_Occupied
    end

    @assert add_vertex!(mg)
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)
    push!(_label2node(lh), label => current_id)

    @assert !is_cyclic(mg)

    return Success
end

function _add_relation!(lh::LabelHierarchy; super::HLabel, sub::HLabel)
    if !_contains(lh, super) || !_contains(lh, sub)
        return Label_Missing
    elseif super == sub
        return Label_Duplication
    elseif _has_relation(lh, super, sub)
        return Relation_Occupied
    end

    mg, super_id, sub_id = _hierarchy(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    @assert add_edge!(mg, sub_id, super_id)

    @assert !is_cyclic(mg)

    return Success
end

function _remove_label!(lh::LabelHierarchy, label::HLabel)
    #@assert label != No_Label
    mg = _hierarchy(lh)
    n = _get_nodeid(lh, label)
    result = rem_vertex!(mg, n)

    result && update_l2n!(lh)
    return result
end

function _remove_relation!(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    mg = _hierarchy(lh)
    n1 = _get_nodeid(lh, label1)
    n2 = _get_nodeid(lh, label2)
    result = rem_edge!(mg, n1, n2) || rem_edge!(mg, n2, n1)

    result && update_l2n!(lh)
    return result
end

function _super_id(lh::LabelHierarchy, id::Integer)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, id::Integer)
    return inneighbors(_hierarchy(lh), id)
end

function _super_id(lh::LabelHierarchy, label::HLabel)
    #@assert label != No_Label
    id = _get_nodeid(lh, label)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, label::HLabel)
    #@assert label != No_Label
    id = _get_nodeid(lh, label)
    return inneighbors(_hierarchy(lh), id)
end

function _super(lh::LabelHierarchy, label::HLabel)
    super_ids = _super_id(lh, label)
    [_get_label(lh, i) for i in super_ids]
end

function _sub(lh::LabelHierarchy, label::HLabel)
    sub_ids = _sub_id(lh, label)
    [_get_label(lh, i) for i in sub_ids]
end

function _issuper(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    _get_nodeid(lh, lhs) ∈ _super_id(lh, rhs)
end

function _issub(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    _get_nodeid(lh, lhs) ∈ _sub_id(lh, rhs)
end

function ==(lhs::LabelHierarchy, rhs::LabelHierarchy)
    # Check Hierarchy(lhs) == Hierarchy(rhs)
    for lhs_label in _labels(lhs)
        if !_contains(rhs, lhs_label)
            return false
        end
        l_super = _super(lhs, lhs_label) |> Set
        r_super = _super(rhs, lhs_label) |> Set
        l_sub   = _super(lhs, lhs_label) |> Set
        r_sub   = _super(rhs, lhs_label) |> Set
        if !(l_super == r_super && l_sub == r_sub)
            return false
        end
    end

    return _label2node(lhs) |> length == _label2node(rhs) |> length
end

function _root_id(lh::LabelHierarchy)
    mg = _hierarchy(lh)
    root_id = filter(i -> isempty(_super_id(lh, i)), 1:nv(mg))

    @assert length(root_id) == 1
    return root_id[1]
end

function _root(lh::LabelHierarchy)
    i = _root_id(lh)
    return _get_label(lh, i)
end

# treeではないのでDFSは使えない
function _depth(lh::LabelHierarchy)
    sub_ids = _sub_id(lh, _root_id(lh))
    depth = 0
    while !isempty(sub_ids)
        depth += 1
        sub_ids = mapreduce(i -> _sub_id(lh, i), append!, sub_ids)
    end

    # excluding atom label
    return depth - 1
end



include("test.jl")
end #module
