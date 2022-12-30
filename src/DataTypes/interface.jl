function heierchy(s::System, hname::Symbol)
    s.heierchy[hname]
end

function heierchies(s::System)
    heierchy(s) |> keys
end

function labels(s::System, hname::Symbol)
    h = heierchy(s, hname)
    HeierchyLabels.labels(h)
end

function add_label!(s::System, hname::Symbol; label::Label, atom_ids::Vector{<:Integer}, super::Label, sub::Label)
    h = heierchy(s, hname)
    HeierchyLabels.add_label!(h, label, atom_ids; super = super, sub = sub)
end

