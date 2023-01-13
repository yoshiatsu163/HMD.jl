module HierarchyLabels

using Graphs
using Match
using MetaGraphs

using ..DataTypes: Id, Category

import Base: getindex
import Base: ==

export Label, LabelHierarchy
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

@enum Result notfound=1

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
function _add_label!(lh::LabelHierarchy, label::Label; super::Label, sub::Label, insert=false)
    mg = _hierarchy(lh)

    _is_add_ready(lh, label; super=super, sub=sub, insert)
    ready_insert = if insert && super != sub != No_Label
        _is_insert_ready(lh, label; super=super, sub=sub)
    elseif insert && (super == No_Label || sub == No_Label)
        error("super or sub is No_Label. ")
    else
        false
    end

    # add node for label
    if !add_vertex!(mg)
        error("failed to add node to Hierarchy graph.")
    end
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)
    push!(_label2node(lh), label => current_id)

    # relation removal
    if ready_insert && !_remove_relation!(lh, sub, super)
        # if _remove_relation! failed
        _remove_label!(lh, label)
        pop!(_label2node(lh), label)
        error("Relation removal between super and sub failed during inserting label. ")
    end
    # add_relation
    if super != No_Label
        _add_relation!(lh; sub=label, super=super)
    end
    if sub != No_Label
        _add_relation!(lh; sub=sub, super=label)
    end
end

function _add_relation!(lh::LabelHierarchy; super::Label, sub::Label)
    if super == No_Label || sub == No_Label
        error("super or sub is No_Label")
    elseif !_contains(lh, super)
        error("super not found ")
    elseif !_contains(lh, sub)
        error("sub not found ")
    elseif super == sub
        error("super and sub are the same ")
    end

    mg, super_id, sub_id = _hierarchy(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    if has_edge(mg, sub_id, super_id)
        error("relation already exists ")
    elseif has_edge(mg, super_id, sub_id)
        error("reversed relation already exists ")
    end

    if !add_edge!(mg, sub_id, super_id)
        error(
            """add relation failed. LabelHierarchy is broken if you called "_add_relation!"
            via "_add_label!" """
        )
    end
end

function _is_add_ready(lh::LabelHierarchy, label::Label; super::Label, sub::Label, insert)
    # Chack if label already exists
    if _contains(lh, label)
        error("label already exists")
    end

    # Check wether super label and sub label exists
    if super != No_Label && !_contains(lh, super)
        error("super not found")
    end
    if sub != No_Label && !_contains(lh, sub)
        error("sub not found")
    end

    if super != sub != No_Label
        # Check if super label and sub label is equal
        if super == sub
            error("super and sub are the same")
        end
        #Check if relation exists
        if !insert && _has_relation(lh, super, sub)
            error(
                """relation between super and sub already exists.
                Set kwarg insert = true to insert label between existing relation. """
            )
        end
    end

    return nothing
end

function _is_insert_ready(lh::LabelHierarchy, label::Label; super::Label, sub::Label)
    if super == No_Label || sub == No_Label
        error("super label and sub label must not be No_Label when insert = true")
    end

    # Check if label has relation between super and sub
    if _contains(lh, label)
        if _has_relation(lh, label, super)
            error("""label and super already has relation""")
        elseif _has_relation(lh, label, sub)
            error("""label and sub already has relation""")
        end
    end

    # Check there is no relation between `super` and `sub`
    ready_insert = if _has_relation(lh, super, sub)
        true
    elseif !_has_relation(lh, super, sub)
        error("relation between super and sub not found.")
    end

    return ready_insert
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
    return rem_edge!(mg, n1, n2)
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
