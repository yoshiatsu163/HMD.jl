Base.@kwdef mutable struct Heierchy
    mg::MetaDiGraph = MetaDiGraph()
    label2node::Dict{Label, Int64} = Dict{Label, Int64}()
end

function _heierchy(h::Heierchy)
    h.mg
end

function _label2node(h::Heierchy)
    h.label2node    
end

function _contains(h::Heierchy, label::Label)
    label ∈ _label2node(h) |> keys
end

function _get_nodeid(h::Heierchy, label::Label)
    _label2node(h)[label]
end

function _get_label(h::Heierchy, id::Integer)
    get_prop(_heierchy(h), id, :label)
end

function _add_label!(h::Heierchy, label::Label; super::Union{Nothing, Label}, sub::Union{Nothing, Label})
    mg = _heierchy(h)

    # chack if label already exists
    if _contains(h, label)
        error("label already exists")
    end

    # check wether super and sub label exists
    if !isnothing(super) && !_contains(h, super)
        error("super not found")
    end
    if !isnothing(sub) && !_contains(h, sub)
        error("sub not found")
    end

    if !isnothing(super) && !isnothing(sub) && super == sub
        error("super and sub are the same")
    end

    # add node for label
    if !add_vertex!(mg)
        error("failed to add node to heierchy graph. ")
    end
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)
    push!(_label2node(h), label => current_id)

    # add_relation
    if !isnothing(super)
        _add_relation!(h; sub = label, super = super)
    end
    if !isnothing(sub)
        _add_relation!(h; sub = sub, super = label)
    end
end

function _add_relation!(h::Heierchy; super::Label, sub::Label)
    if !_contains(h, super)
        error("super not found ")
    elseif !_contains(h, sub)
        error("sub not found ")
    elseif super == sub
        error("super and sub are the same ")
    end

    mg, super_id, sub_id = _heierchy(h), _get_nodeid(h, super), _get_nodeid(h, sub)
    if has_edge(mg, sub_id, super_id)
        error("relation already exists ")
    elseif has_edge(mg, super_id, sub_id)
        error("reversed relation already exists ")
    end

    if !add_edge!(mg, sub_id, super_id)
        error("add relation failed")
    end
end

function _super(h::Heierchy, label::Label)
    id = _get_nodeid(h, label)
    super_ids = outneighbors(_heierchy(h), id)
    [_get_label(h, i) for i in super_ids]
end

function _sub(h::Heierchy, label::Label)
    id = _get_nodeid(h, label)
    sub_ids = inneighbors(_heierchy(h), id)
    [_get_label(h, i) for i in sub_ids]
end

function _issuper(h::Heierchy, lhs::Label, rhs::Label)
    lhs ∈ _super(h, rhs)
end

function _issub(h::Heierchy, lhs::Label, rhs::Label)
    lhs ∈ _sub(h, rhs)
end
