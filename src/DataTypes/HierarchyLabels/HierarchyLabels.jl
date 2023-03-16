module HierarchyLabels

using Graphs
using MLStyle

#using ..DataTypes: Id, Category, StaticString, name, SerializedCategory
import ..DataTypes: serialize, deserialize

import Base: getindex
import Base: ==
import Base: ∈
import Base: ∋
import Base: ∉
import Base: show

export HLabel, H_Label, LabelHierarchy
export LabelResult, Label_Missing, Label_Occupied, Label_Duplication, Relation_Missing, Relation_Occupied, Success
export id, type, ==
export _add_label!, _add_labels!, _add_relation!, _remove_label!, _remove_relation!
export _label2node, _contains, ∈, ∋, ∉, _has_relation ,_get_nodeid, getindex
export _issuper, _issub, _super_id, _sub_id, _super, _sub
export _root, _depth
export PackedHierarchy, serialize, deserialize

# AbstractLabel must define constructor s.t. ALabel(id::Integer, type)
abstract type AbstractLabel end

struct HLabel <: AbstractLabel
    type::String
    id::Int64
end

function id(label::HLabel)
    return label.id
end

function type(label::HLabel)
    return label.type
end

#@enum LabelResult Label_Missing=1 Label_Occupied=2 Label_Duplication=3 Relation_Missing=4 Relation_Occupied=5 Success=6
@data LabelResult begin
    Label_Missing
    Label_Occupied
    Label_Duplication
    Relation_Missing
    Relation_Occupied
    Success
end

function ==(lhs::HLabel, rhs::HLabel)
    id(lhs) == id(rhs) && type(lhs) == type(rhs)
end

function Base.show(io::IO, ::MIME"text/plain", label::HLabel)
    i = id(label)
    t = type(label)
    print(io, "HLabel(\"$t\", $i)")
end

function Base.show(label::HLabel)
    i = id(label)
    t = type(label)
    print("HLabel(\"$t\", $i)")
end

function Base.print_to_string(label::HLabel)
    i = id(label)
    t = type(label)
    return "HLabel(\"$t\", $i)"
end

Base.@kwdef mutable struct LabelHierarchy
    g::DiGraph{Int64} = DiGraph{Int64}()
    labels::Vector{HLabel} = Vector{HLabel}()
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
    return lh.g
end

function _label2node(lh::LabelHierarchy)
    return lh.label2node
end

function _labels(lh::LabelHierarchy)
    return lh.labels
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
    return _labels(lh)[id]
end

function _set_label!(lh::LabelHierarchy, label::HLabel, id::Integer)
    _labels(lh)[id] = label
    _label2node(lh)[label] = id
end

function _add_label!(lh::LabelHierarchy, label::HLabel)
    g = _hierarchy(lh)

    if _contains(lh, label)
        return Label_Occupied
    end

    @assert add_vertex!(g)
    current_id = nv(g)
    push!(_labels(lh), label)
    merge!(_label2node(lh), Dict(label => current_id))

    @assert !is_cyclic(g)

    return Success
end

function _add_labels!(lh::LabelHierarchy, labels::AbstractVector{HLabel})
    g = _hierarchy(lh)

    new_range = (nv(g) + one(nv(g))) : (nv(g) + length(labels))
    @assert add_vertices!(g, length(labels)) == length(labels)
    append!(lh.labels, labels)
    merge!(_label2node(lh), Dict(label => node_id for (label, node_id) in zip(labels, new_range)))

    return Success
end

function _add_relation!(lh::LabelHierarchy; super::HLabel, sub::HLabel, unsafe::Bool=false)
    if !unsafe
        if !_contains(lh, super) || !_contains(lh, sub)
            return Label_Missing
        elseif super == sub
            return Label_Duplication
        elseif _has_relation(lh, super, sub)
            return Relation_Occupied
        end
    end

    g, super_id, sub_id = _hierarchy(lh), _get_nodeid(lh, super), _get_nodeid(lh, sub)
    @assert add_edge!(g, sub_id, super_id)

    if !unsafe
        @assert !is_cyclic(g)
    end

    return Success
end

function _remove_label!(lh::LabelHierarchy, label::HLabel)
    if !_contains(lh, label)
        error("LabelHierarchy does not contain $label. ")
    end
    g, id = _hierarchy(lh), _get_nodeid(lh, label)

    @assert rem_vertex!(g, id)
    # when removing a vertex, the vertex id of the last vertex id is changed to keep vertex ids consecutive.
    end_label = pop!(_labels(lh))
    delete!(_label2node(lh), label)
    _set_label!(lh, end_label, id)

    return nothing
end

function _remove_relation!(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !_contains(lh, label1) || !_contains(lh, label2)
        error("LabelHierarchy does not contain $label1 or $label2. ")
    end

    g = _hierarchy(lh)
    n1, n2 = _get_nodeid(lh, label1), _get_nodeid(lh, label2)

    @assert has_edge(g, n2, n1)
    @assert rem_edge!(g, n1, n2) || rem_edge!(g, n2, n1)

    return nothing
end

function _contains(lh::LabelHierarchy, label::HLabel)
    return haskey(_label2node(lh), label)
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

function _has_relation(lh::LabelHierarchy, label1::HLabel, label2::HLabel)
    if !(_contains(lh, label1) && _contains(lh, label2))
        error("labels not found in LabelHierarchy. ")
    end
    return _issuper(lh, label1, label2) || _issub(lh, label1, label2)
end

function _super_id(lh::LabelHierarchy, id::Integer)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, id::Integer)
    return inneighbors(_hierarchy(lh), id)
end

function _super_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return outneighbors(_hierarchy(lh), id)
end

function _sub_id(lh::LabelHierarchy, label::HLabel)
    id = _get_nodeid(lh, label)
    return inneighbors(_hierarchy(lh), id)
end

function _super(lh::LabelHierarchy, label::HLabel)
    super_ids = _super_id(lh, label)
    if isempty(super_ids)
        return Vector{HLabel}(undef, 0)
    else
        return view(_labels(lh), super_ids)
    end
end

function _sub(lh::LabelHierarchy, label::HLabel)
    sub_ids = _sub_id(lh, label)
    if isempty(sub_ids)
        return Vector{HLabel}(undef, 0)
    else
        return view(_labels(lh), sub_ids)
    end
end

function _issuper(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) ∈ _super_id(lh, rhs)
end

function _issub(lh::LabelHierarchy, lhs::HLabel, rhs::HLabel)
    return _get_nodeid(lh, lhs) ∈ _sub_id(lh, rhs)
end

function ==(lhs::LabelHierarchy, rhs::LabelHierarchy)
    if _label2node(lhs) != _label2node(rhs)
        return false
    end

    if Set(_labels(lhs)) != Set(_labels(rhs))
        return false
    end

    lg, rg = _hierarchy(lhs), _hierarchy(rhs)
    for lhs_id in vertices(lg)
        lhs_label = _get_label(lhs, lhs_id)
        super_labels = _super(lhs, lhs_label) |> Set
        sub_labels = _sub(lhs, lhs_label) |> Set
        rhs_label = _get_label(rhs, _get_nodeid(rhs, lhs_label))
        if super_labels != Set(_super(rhs, rhs_label)) || sub_labels != Set(_sub(rhs, rhs_label))
            return false
        end
    end

    return true
end

function _root_id(lh::LabelHierarchy)
    g = _hierarchy(lh)
    root_id = filter(i -> isempty(_super_id(lh, i)), 1:nv(g))

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

struct PackedHierarchy
    num_node::Int64
    edges_org::Vector{Int64}
    edges_dst::Vector{Int64}
    label_ids::Vector{Int64}
    # splitting Vector{String}
    chars::Vector{UInt8}
    bounds::Vector{Int64}
end

function serialize(lh::LabelHierarchy)
    # split Vecotr{HLabel} to Vector{Int64} and Vector{String}
    label_ids = [id(label) for label in _labels(lh)]
    chars, bounds = serialize([type(label) for label in _labels(lh)])

    g = _hierarchy(lh)
    num_node = nv(g)
    edges_org, edges_dst = begin
        orig, dest = Vector{Int64}(undef, ne(g)), Vector{Int64}(undef, ne(g))
        for (i, edge) in enumerate(edges(g))
            orig[i], dest[i] = src(edge), dst(edge)
        end
        orig, dest
    end

    return PackedHierarchy(num_node, edges_org, edges_dst, label_ids, chars, bounds)
end

function deserialize(ph::PackedHierarchy)
    g = begin
        g = DiGraph()
        add_vertices!(g, ph.num_node)
        for (o, d) in zip(ph.edges_org, ph.edges_dst)
            add_edge!(g, o, d)
        end
        g
    end

    label_ids = ph.label_ids
    label_types = deserialize(ph.chars, ph.bounds)
    labels = [HLabel(type, id) for (type, id) in zip(label_types, label_ids)]

    lh = LabelHierarchy()
    lh.g = g
    lh.labels = labels
    lh.label2node = Dict(labels .=> 1:length(labels))

    return lh
end

include("test.jl")

end #module
