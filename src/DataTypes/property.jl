function prop_names(s::System)
    s.props |> keys
end

function props(s::System, pname::AbstractString)
    s.props[pname]
end

function prop(s::System, pname::AbstractString, label::Label)
    props(s, pname)[label]
end

function labels_in_prop(s::System, pname::AbstractString)
    props(s, pname) |> keys
end

function ∈(label::Label, prp::Dict{Label, Any})
    label ∈ keys(prp)
end

function ∋(prp::Dict{Label, Any}, label::Label)
    label ∈ keys(prp)
end

function ∉(label::Label, prp::Dict{Label, Any})
    label ∉ keys(prp)
end

function add_prop!(s::System, pname::AbstractString)
    push!(s.props, pname => Dict{Label, Any}())
end

function add_prop!(s::System, pname::AbstractString, label::Label, p::Any)
    prp = props(s, pname)
    if label ∈ prp
        error("label $(label) already exists. ")
    end
    push!(prp, label => p)
end

function set_prop!(s::System, pname::AbstractString, label::Label, p::Any)
    if label ∉ props(s, pname)
        error("label $(label) not found in property $(pname). ")
    end
    props(s, pname)[label] = p
end
