module HeierchyLabels

using MetaGraphs

export Label, LabelHeierchy, Tree
export super, sub

struct Label
    name::String
end

const Tree{T} = MetaDiGraph{T} where T <: Integer

mutable struct LabelHeierchy{T <: Integer, S <: AbstractString}
    labels::Tree{T}
    graphs::Dict{Label, MetaGraph}
    name::S
end

function root(t::Tree)
    !is_connected && error("Tree is disconnetted. ")
    !is_cyclic && error("Tree has cyclic parts. ")

    climb(g, v) = begin
        p = inneighbors(t, v)
        if length(p) == 1
            climb(g, v)
        elseif length(p) == 0
            return v
        else
            length(p) >= 2 && error("Node $v has multiple parents. ")
        end
    end

    climb(t.tree, 1)
end

function find_node(t::Tree, l::Label)
    filter(v->props(v)[:label]==l, vertices(t))
end

function super(t::Tree, l::Label)
    list = map(v->inneighbors(t, v), find_node(t, l))
    unique(list)
end

function sub(t::Tree, l::Label)
    list = map(v->outneighbors(t, v), find_node(t, l))
    unique(list)
end

function rename!(t::Tree, l::Label, new::AbstractString)
    v = find_node(t, l)
    prop = deepcopy(props(v))
    prop[:label].name = new
    set_props!(t, v, prop)
end

#module
end