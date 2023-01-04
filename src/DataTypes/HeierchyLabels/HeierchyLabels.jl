module HeierchyLabels

using Graphs, MetaGraphs
import Base: getindex

export Label, LabelHeierchy
export id, type, l2a, add_atom!, add_label!, add_relation!, super, sub, issuper, issub, getindex
export Label2atom
export NoLabel

using ..DataTypes: Id

struct Label
    id::Id{Label}
    type::String
end

const NoLabel = Label(typemin(Int64), "NoLabel")

function id(label::Label)
    label.id
end

function type(label::Label)
    label.type
end

include("label2atom.jl")
include("heierchy.jl")

Base.@kwdef mutable struct LabelHeierchy
    heierchy::Heierchy = Heierchy()
    l2a::Label2atom = Label2atom()
end

function heierchy(lh::LabelHeierchy)
    lh.heierchy
end

function l2a(lh::LabelHeierchy)
    lh.l2a
end

function add_atom!(lh::LabelHeierchy, n::Integer)
    _add_atom!(l2a(lh), n)
end

function add_label!(lh::LabelHeierchy, label::Label, atom_ids::Vector{<:Integer}; super::Label, sub::Label)
    _add_label!(heierchy(lh), label; super = super, sub = sub)
    _add_label!(l2a(lh), label, atom_ids)
end

function add_relation!(lh::LabelHeierchy; super::Label, sub::Label)
    h = heierchy(lh)
    _add_relation!(h, super = super, sub = sub)
end

function super(lh::LabelHeierchy, label::Label)
    h = heierchy(lh)
    _super(h, label), l2a(lh)[label]
end

function sub(lh::LabelHeierchy, label::Label)
    h = heierchy(lh)
    _sub(h, label), l2a(lh)[label]
end

function issuper(lh::LabelHeierchy, lhs::Label, rhs::Label)
    rhs ∈ super(lh, lhs)[1]
end

function issub(lh::LabelHeierchy, lhs::Label, rhs::Label)
    rhs ∈ sub(lh, lhs)[1]
end

function labels(lh::LabelHeierchy)
    l2a(lh) |> _labels
end

function getindex(lh::LabelHeierchy, label::Label)
    l2a(lh)[label]
end

function getindex(lh::LabelHeierchy, atomid::Integer)
    l2a(lh)[atomid]
end

include("test.jl")
end #module