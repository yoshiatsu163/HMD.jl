module lammps

using Pipe

@enum LammpsStyle other=0 atomic=1 bond=2 angle=3 full=4

mutable struct Box
    xdim::String
    ydim::String
    zdim::String
    triclinic::String
end

mutable struct Header
    comment::String
    n_attributes::Vector{Int64}
    n_types::Vector{Int64}
    box::Box
end

mutable struct Sections
    header::Vector{String}
    topologies::Vector{String}
    atoms::Vector{String}
    velocities::Vector{String}
    bonds::Vector{String}
    angles::Vector{String}
    diheds::Vector{String}
    imprs::Vector{String}

    Sections() = new()
end

# MetaGraphへの変換を優先してangle以降はいったん無視
# labelの仕様がうまく定義できたらangle以降も読み取る
# よって以下はいったん破棄
# lines_between(signal1, signal2)を定義すると楽かも
function readfile(filename::String)
    fp=open(filename, "r")
    lines = readlines(fp)
    close(fp)
    # check file
    if !("Atoms" in lines)
        error("There must be section \"Atoms\" in lammps data file.")
    elseif !("Masses" in lines)
        error("There must be section \"Masses\" in lammps data file.")
    end

    # check lammps style
    lammps_style = other
    if "Atoms" in lines
        lammps_style = atomic
    elseif "Atoms" in lines && "Bonds" in lines
        lammps_style = bond
    elseif "Atoms" in lines && "Bonds" in lines && "Angles" in lines
        lammps_style = angle
        println("WARNING: reading style \"angle\" as style \"bond\".")
    elseif "Atoms" in lines && "Bonds" in lines && "Angles" in lines && "Dihedrals" in lines
        lammps_style = full
        println("WARNING: reading style \"full\" as style \"bond\".")
    end
    println("Detected atom style: $lammps_style")
    if lammps_style == atomic
        error("Atomic style is not currently supported.")
    elseif "Velocities" in lines
        println("WARNING: Section \"Velocities\" is currently ignored.")
    end

    section = Sections()
    # header section
    header_end  = findfirst(l->occursin("Masses", l), lines)
    save2section!(section, :header, lines[1:header_end-1])

    # atom section
    atoms_top   = findfirst(l->occursin("Atoms", l), lines)
    atoms_end   = findfirst(l->occursin(r"Velocities|Bonds", l), lines)

    # bond section
    bonds_top   = findfirst(l->occursin("Bnods", l), lines)
    bonds_end   = findfirst(l->occursin("Angles", l), lines)

    save2section!(section, :header, lines[atoms_top+1:atoms_end-1])
    save2section!(section, :header, lines[bonds_top+1:bonds_end-1])
    
end

function readfile(filename::String)
    fp=open(filename, "r")
    lines = readlines(fp)
    close(fp)
    # check file
    if !("Atoms" in lines)
        error("There must be section \"Atoms\" in lammps data file.")
    elseif !("Masses" in lines)
        error("There must be section \"Masses\" in lammps data file.")
    end

    # check lammps style
    lammps_style = other
    if "Atoms" in lines
        lammps_style = atomic
    elseif "Atoms" in lines && "Bonds" in lines
        lammps_style = bond
    elseif "Atoms" in lines && "Bonds" in lines && "Angles" in lines
        lammps_style = angle
    elseif "Atoms" in lines && "Bonds" in lines && "Angles" in lines && "Dihedrals" in lines
        lammps_style = full
    end
    println("Detected atom style: $lammps_style")
    
    section = Sections()
    # header section
    header_end  = findfirst(l->occursin("Masses", l), lines)
    save2section!(section, :header, lines[1:header_end-1])

    # atom section
    atoms_top   = findfirst(l->occursin("Atoms", l), lines)
    atoms_end   = findfirst(l->occursin(r"Velocities|Bonds", l), lines)

    # bond section
    bonds_top   = findfirst(l->occursin("Bnods", l), lines)
    bonds_end   = findfirst(l->occursin("Angles", l), lines)

    # angle section
    angles_top  = bonds_end
    angles_end  = findfirst(l->occursin("Dihedrals", l), lines)

    # dihedral section
    diheds_top  = angles_end
    diheds_end  = findfirst(l->occursin("Impropers", l), lines)

    save2section!(section, :atoms, lines[atoms_top+1:atoms_end-1])
    if lammps_style == bond
        save2section!(section, :bonds, lines[bonds_top+1:end])
    elseif lammps_style == angle
        save2section!(section, :bonds, lines[bonds_top+1:bonds_end-1])
        save2section!(section, :angles, lines[angles_top+1:end])
    elseif lammps_style == full
        save2section!(section, :bonds, lines[bonds_top+1:bonds_end-1])
        save2section!(section, :angles, lines[angles_top+1:angles_end-1])
        if diheds_end == nothing
            save2section!(section, :diheds, lines[angles_top+1:end])
        else
            save2section!(section, :diheds, lines[diheds_top+1:diheds_end-1])
            save2section!(section, :imprs, lines[diheds_end+1:end])
        end
    end

    sections
end

function write_lmpdat(filename, atoms, topologies, coeffs, box)
    fp = open(filename, "w")

    #header
    write(fp, "# This data file is created by manyMolecule.jl\n\n")
    @pipe select(atoms, :atomid) |> maximum |> write(fp, "$_ atoms\n")
    @pipe filter(t->t[:category]=="Bonds"    , topologies) |> length(select(_, :tid)) |> write(fp, "$_ bonds\n")
    @pipe filter(t->t[:category]=="Angles"   , topologies) |> length(select(_, :tid)) |> write(fp, "$_ angles\n")
    @pipe filter(t->t[:category]=="Dihedrals", topologies) |> length(select(_, :tid)) |> write(fp, "$_ dihedrals\n")
    @pipe filter(t->t[:category]=="Impropers", topologies) |> length(select(_, :tid)) |> write(fp, "$_ impropers\n")
    write(fp, "\n")

    @pipe select(atoms, :atomtype) |> maximum |> write(fp, "$_ atom types\n")
    @pipe filter(t->t[:category]=="Bond Coeffs"    , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ bond types\n")
    @pipe filter(t->t[:category]=="Angle Coeffs"   , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ angle types\n")
    @pipe filter(t->t[:category]=="Dihedral Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ dihedral types\n")
    @pipe filter(t->t[:category]=="Improper Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ improper types\n")
    write(fp, "\n")
    
    write(fp, "0.0 $(box[1]) xlo xhi\n")
    write(fp, "0.0 $(box[2]) ylo yhi\n")
    write(fp, "0.0 $(box[3]) zlo zhi\n")
    write(fp, "0.0 0.0 0.0 xy xz yz\n")

    write(fp, "\nMasses\n\n")
    for molname in select(atoms, :molname) |> unique
        for row in filter(t->(t[:category]=="Masses" && t[:molname]==molname), coeffs)
            write(fp, "$(row[:cid])  $(row[:coeff])\n")
        end
    end
    write(fp, "\n")

    for cat in @pipe select(coeffs, :category) |> unique |> filter(!=("Masses"),_)
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, coeffs)
            write(fp, "$(row[:cid]) $(row[:style]) $(row[:coeff])\n")
        end
        write(fp, "\n")
    end

    write(fp, "Atoms\n\n")
    for row in atoms
        write(fp, "$(row[:atomid])  $(row[:molid])  $(row[:atomtype])    $(row[:q])     $(row[:x][1])     $(row[:x][2])     $(row[:x][3])\n")
    end
    write(fp, "\n")

    for cat in select(topologies, :category) |> unique
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, topologies)
            str = "$(row[:tid]) $(row[:cid])"
            for aid in row[:atomid]
                str *= "     $(aid)"
            end
            write(fp, str * "\n")
        end
        write(fp, "\n")
    end

    close(fp)
end

function save2section!(section, field, lines)
    @pipe filter(!isempty, lines[s:f-1]) |> map(l->replace(l, "\t"=>"    "), _) |> setfield!(sections, field, _)
end

function in(needle::Union{AbstractString,AbstractPattern,AbstractChar}, haystack::Vector{String})
    findfirst(l->occursin(needle, l), haystack) != nothing
end

end