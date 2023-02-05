function prop_names(s::AbstractSystem)
    s.props |> keys
end

function props(s::AbstractSystem, pname::AbstractString)
    s.props[pname]
end

function prop(s::AbstractSystem, pname::AbstractString, label::HLabel)
    props(s, pname)[label]
end

function labels_in_prop(s::AbstractSystem, pname::AbstractString)
    props(s, pname) |> keys
end

function ∈(label::HLabel, prp::Dict{HLabel, Any})
    label ∈ keys(prp)
end

function ∋(prp::Dict{HLabel, Any}, label::HLabel)
    label ∈ keys(prp)
end

function ∉(label::HLabel, prp::Dict{HLabel, Any})
    label ∉ keys(prp)
end

function add_prop!(s::AbstractSystem, pname::AbstractString)
    push!(s.props, pname => Dict{HLabel, Any}())
end

function add_prop!(s::AbstractSystem, pname::AbstractString, label::HLabel, p::Any)
    prp = props(s, pname)
    if label ∈ prp
        error("label $(label) already exists. ")
    end
    push!(prp, label => p)
end

function set_prop!(s::AbstractSystem, pname::AbstractString, label::HLabel, p::Any)
    if label ∉ props(s, pname)
        error("label $(label) not found in property $(pname). ")
    end
    props(s, pname)[label] = p
end
