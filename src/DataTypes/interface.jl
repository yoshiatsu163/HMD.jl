function heiherchy()
    
end

function set_heiherchy!()
    
end

# labelname => atomidに限るかlabelname => labelも可能にするか
function set_labels!(s::System, name::Symbol, map::AbstractDict; reverse = false)
    map_name = (name, :aid)
    mm = labels(s)
    set_map!(mm, map_name, map)
    if reverse
        add_reverse!(mm, ())
    end
end