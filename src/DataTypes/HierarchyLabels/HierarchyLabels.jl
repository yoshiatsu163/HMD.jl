module HierarchyLabels

using Graphs
using Match
using MetaGraphs

using ..DataTypes: Id, Category

import Base: getindex
import Base: ==

export Label, LabelHierarchy, LabelResult
export id, type, ==
export _add_label!, _add_relation!, _remove_label!, _remove_relation
export _label2node, _contains, _has_relation ,_get_nodeid, getindex
export _issuper, _issub, _super_id, _sub_id, _super, _sub

struct Label
    id::Id{Label}
    type::Category{Label}
end

function Label(id::Integer, type::AbstractString)
    Label(Id{Label}(id), Category{Label}(type))
end

const No_Label = Label(typemin(Int64), "No_Label")

@enum LabelResult Label_Missing=1 Label_Occupied=2 Label_Duplication=3 Relation_Missing=4 Relation_Occupied=5 Success=6

function id(label::Label)
    label.id
end

function type(label::Label)
    label.type
end

function ==(lhs::Label, rhs::Label)
    id(lhs) == id(rhs) && type(lhs) == type(rhs)
end

Base.@kwdef mutable struct LabelHierarchy
    mg::MetaDiGraph = MetaDiGraph()
    label2node::Dict{Label, Int64} = Dict{Label, Int64}()
end

function _hierarchy(lh::LabelHierarchy)
    return lh.mg
end

function _label2node(lh::LabelHierarchy)
    return lh.label2node
end

function _contains(lh::LabelHierarchy, label::Label)
    return label ∈ _label2node(lh) |> keys
end

function ∈(label::Label, lh::LabelHierarchy)
    return _contains(lh, label)
end

function ∋(lh::LabelHierarchy, label::Label)
    return _contains(lh, label)
end

function ∉(label::Label, lh::LabelHierarchy)
    return !_contains(lh, label)
end

function _get_nodeid(lh::LabelHierarchy, label::Label)
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
function _has_relation(lh::LabelHierarchy, label1::Label, label2::Label)
    @assert label1 != label2 != No_Label
    if !(_contains(lh, label1) && _contains(lh, label2))
        error("labels not found in LabelHierarchy. ")
    end
    return _issuper(lh, label1, label2) || _issub(lh, label1, label2)
end

# ここでinsertionを定義するのは複雑なのでdatatypesのほうがよさそう
function _add_label!(lh::LabelHierarchy, label::Label)
    mg = _hierarchy(lh)

    if _contains(lh, label)
        return Label_Occupied
    end

    @assert !add_vertex!(mg)
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)
    push!(_label2node(lh), label => current_id)

    return Success
end

function _add_relation!(lh::LabelHierarchy; super::Label, sub::Label)
    if !_contains(lh, super) || !_contains(lh, sub)
        return Label_Missing
    elseif super == sub
        return Label_Duplication
    elseif _has_relation(lh, super, sub)
        return Relation_Occupied
    end

    mg, super_id, sub_id = _hierarchy(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    @assert !add_edge!(mg, sub_id, super_id)

    return Success
end

function _remove_label!(lh::LabelHierarchy, label::Label)
    @assert label != No_Label
    mg = _hierarchy(lh)
    n = _get_nodeid(lh, label)
    return rem_vertex!(mg, n)
end

function _remove_relation!(lh::LabelHierarchy, label1::Label, label2::Label)
    @assert label1 != label2 != No_Label
    mg = _hierarchy(lh)
    n1 = _get_nodeid(lh, label1)
    n2 = _get_nodeid(lh, label2)
    return rem_edge!(mg, n1, n2) || rem_edge!(mg, n2, n1)
end

function _super_id(lh::LabelHierarchy, label::Label)
    @assert label != No_Label
    id = _get_nodeid(lh, label)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, label::Label)
    @assert label != No_Label
    id = _get_nodeid(lh, label)
    return inneighbors(_hierarchy(lh), id)
end

function _super(lh::LabelHierarchy, label::Label)
    super_ids = _super_id(lh, label)
    [_get_label(lh, i) for i in super_ids]
end

function _sub(lh::LabelHierarchy, label::Label)
    sub_ids = _sub_id(lh, label)
    [_get_label(lh, i) for i in sub_ids]
end

function _issuper(lh::LabelHierarchy, lhs::Label, rhs::Label)
    _get_nodeid(lh, lhs) ∈ _super_id(lh, rhs)
end

function _issub(lh::LabelHierarchy, lhs::Label, rhs::Label)
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

include("test.jl")
end #module
