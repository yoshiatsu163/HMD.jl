#TODO: _add_label 周りのロジックが複雑 & エラー発生時のデータ不変性が未確保
#TODO: _add_label ボイラープレートも多い
#TODO: _add_label genericなエラーチェック関数を作る
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

function _has_relation(h::Heierchy, label1::Label, label2::Label)
    issuper(h, label1, label2) || issub(h, label1, label2)
end

function _add_label!(h::Heierchy, label::Label; super::Label, sub::Label, insert=false)
    mg = _heierchy(h)

    # chack if label already exists
    if _contains(h, label)
        error("label already exists")
    end

    # check wether super and sub label exists
    if super != No_Label && !_contains(h, super)
        error("super not found")
    end
    if sub != No_Label && !_contains(h, sub)
        error("sub not found")
    end

    if super != No_Label && sub != No_Label && super == sub
        error("super and sub are the same")
    end

    # Check there is no relation between `super` and `sub`
    ready_insert = if insert && _has_relation(h, super, sub)
        true
    elseif insert && !_has_relation(h, super, sub)
        error("relation between super and sub not found.")
    elseif !insert && _has_relation(h, super, sub)
        error(
            """relation between super and sub alreeady exists.
            Set kwarg insert = true to insert label between existing relation. """
        )
    else
        false
    end

    # add node for label
    if !add_vertex!(mg)
        error("failed to add node to heierchy graph.")
    end
    current_id = nv(mg)
    set_prop!(mg, current_id, :label, label)
    push!(_label2node(h), label => current_id)

    # add_relation
    if ready_insert && !_remove_relation!(h, super, sub)
        _remove_label!(h, label)
        pop!(_label2node(h), label)
        error("Relation removal between super and sub failed. ")
    end
    if super != No_Label
        _add_relation!(h; sub=label, super=super)
    end
    if sub != No_Label
        _add_relation!(h; sub=sub, super=label)
    end
end

function _add_relation!(h::Heierchy; super::Label, sub::Label)
    if super == No_Label || sub == No_Label
        error("super or sub is No_Label")
    elseif !_contains(h, super)
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

function _remove_label!(h::Heierchy, label::Label)
    mg = _heierchy(h)
    n = _get_nodeid(h, label)
    return rem_vertex!(mg, n)
end

function _remove_relation!(h::Heierchy, label1::Label, label2::Label)
    mg = _heierchy(h)
    n1 = _get_nodeid(h, label1)
    n2 = _get_nodeid(h, label2)
    return rem_edge!(mg, n1, n2)
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

function ==(lhs::Heierchy, rhs::Heierchy)
    _heierchy(lhs) == _heierchy(rhs) && _label2node(lhs) == _label2node(rhs)
end
