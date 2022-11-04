module DataStructure

 
export BoundingBox, Topology, TimeSeries, Trajectory

using StaticArrays, Graphs, MetaGraphs

include("HeierchyLabels.jl")
using .HeierchyLabels


"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""
struct BoundingBox{T <: AbstractFloat}
    xvec::SVector{3, T}
    yvec::SVector{3, T}
    zvec::SVector{3, T}
end

macro pos(ex)
    ex.head != :(.) && error("""Expression must be like "arr[i].x" """)
    ex.args[1].head != :(ref) && error!("""Expression must be like "arr[i].x" """)
    
    array_name = ex.args[1].args[1]
    col   = ex.args[1].args[2]
    row = if ex.args[2] == :(:x)
        :(1)
    elseif ex.args[2] == :(:y)
        :(2)
    elseif ex.args[2] == :(:z)
        :(3)
    else
        error("""Symbol after "." must be "x" or "y" or "z" """)
    end

    Expr(:ref, array_name, row, col)
end

#const Position = SVector{3, T <: AbstractVector}

mutable struct TimeSeries{T <: Real}
    time::Vector{T}
    prop::Vector
end

mutable struct Topology{T <: AbstractFloat, U <: AbstractRange{T}}
    graph::Vector{MetaGraph}
    #labels::Dict{String, Vector{LabelHeierchy}}
    #lifetime::Vector{U}
end

function ∈(lhs::Real, rhs::AbstractRange)
    first(rhs) <= lhs <= last(rhs)
end

function (topo::Topology)(time)
    ind = findfirst(l -> time∈l, topo.lifetime)
    @assert !isnothing(ind)
    Topology(topology.graph[ind], topology.labels[ind], time:time)
end

mutable struct Trajectory{T <: Integer, U <: AbstractFloat, X <: Any}
    time::Vector{U}
    topology::Topology{U} #with atomid, charge, element, LabelHeierchy 1本になるはず
    prop::Dict{Symbol, X} # labelの変化にはどう対応するか？
    position::Vector{Matrix{U}}
end

mutable struct System
    time
    topology
    prop
    position
end

mutable struct Traj
    time::Vector{Int}
    system::Vector{System}
    prop::Dict{Symbol, X}
end

#function Trajectory(g::MetaGraph)
#    if !all(i -> keys(props(g, i)) ∋ :position , vertices(g))
#        error("Each vertices in the Metagraph must have :positon property")
#    end
#    position = mapreduce(i -> props(g, i)[:position], hcat, vertices(g))
#    Trajectory( [1.0],
#                Topology(g, Tree(), 0.0:0.0),
#                Dict(),
#                position)
#end


function select(traj, symbol; step, mol, atom, label_condition)
    #filter(x, view(traj, step))
end

"""
    heirechy labelの書き込み
"""
function recieve(traj)
    
end


end #module