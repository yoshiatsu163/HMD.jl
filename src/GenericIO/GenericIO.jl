module GenericIO
    export read
    include("lammps.jl")
    include("mol.jl")
    include("HMDbinary.jl")
    include("HMDtext.jl")

    function read(filename, filetype)
        filename = replace(filename, "\\"=>"/")
        #g = MetaGraph()
        Meta.parse("g = $(filetype).readfile( \"$(filename)\" )") |> eval

        return g
    end
end