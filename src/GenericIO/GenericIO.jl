module GenericIO

export readfile, readfile_MG

include("util.jl")
using .util

include("lammps.jl")
include("mol.jl")
include("HMDbinary.jl")
include("HMDtext.jl")

using ..DataTypes

function readfile_MG(filename, filetype)
    filename = replace(filename, "\\"=>"/")
    #g = MetaGraph()
    Meta.parse("g = $(filetype).readfile( \"$(filename)\" )") |> eval

    return g
end

function readfile(filename, filetype)
    readfile_MG(filename, filetype) |> System
end

end #module