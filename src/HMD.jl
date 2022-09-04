module HMD

using PeriodicTable, StaticArrays, Pipe

include("mol.jl")

function read(filename, filetype)
    $(filetype).readfile( $(filename) )
end

end
