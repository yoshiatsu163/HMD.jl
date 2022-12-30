Base.@kwdef mutable struct Heierchy
    mg::MetaDiGraph = MetaDiGraph()
end

function _heierchy(h::Heierchy)
    h.mg
end

function _get_nodeid(h::Heierchy, label::Label; allow_NA = true)
    mg = _heierchy(h)
    ids = findall(id -> props(mg, id)[:label] == label, vertices(mg))
    if !allow_NA && isempty(ids)
        error("""Label "$(type(label))" not found in the heierchy """)
    else
        return ids
    end
end

function _add_label!(h::Heierchy, label::Label; super::Label, sub::Label)
    mg = _heierchy(h)

    # add node for label
    if !add_vertex!(mg)
        error("failed to add node to heierchy graph. ")
    end
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)

    if super != "nothing"
        super_id   = _get_nodeid(h, super; allow_NA = false)
        @assert length(current_id) == length(super_id) == 1
        add_edge!(h, current_id, super_id)
    end
    if sub != "nothing"
        sub_id     = _get_nodeid(lh, sub; allow_NA = false)
        @assert length(current_id) == length(sub_id) == 1
        add_edge!(h, sub_id, current_id)
    end   
end

function _add_relation!(h::Heierchy, label::String; super::Label, sub::Label)
    if super == sub == "nothing"
        error("""Both super and sub are "nothing" """)
    end

    mg = _heierchy(h)
    if super != "nothing"
        current_id = _get_nodeid(lh, label; allow_NA = false)
        super_id   = _get_nodeid(lh, super; allow_NA = false)
        @assert length(current_id) == length(super_id) == 1
        add_edge!(h, current_id, super_id)
    end
    if sub != "nothing"
        current_id = _get_nodeid(lh, label; allow_NA = false)
        sub_id     = _get_nodeid(lh, sub; allow_NA = false)
        @assert length(current_id) == length(sub_id) == 1
        add_edge!(h, sub_id, current_id)
    end
end

function _super(h::Heierchy, label::Label)
    id = _get_nodeid(h, label)
    outneighbors(_heierchy(h), id)
end

function _sub(h::Heierchy, label::Label)
    id = _get_nodeid(h, label)
    inneighbors(_heierchy(h), id)
end

function _issuper(h::Heierchy, lhs::Label, rhs::Label)
    rhs ∈ _super(h, lhs)
end

function _issub(h::Heierchy, lhs::Label, rhs::Label)
    rhs ∈ _sub(h, lhs)
end
