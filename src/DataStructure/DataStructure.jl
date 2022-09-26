module DataStrucutre

"""
box読み書き
HMDAnalyzeからは内部構造が見えないので一般に読み書きができるapiをここで整備する
"""

export Trajectory, System, Molecule, AbstractStructure, AbstractProperty, Property, select
export Box

abstract type AbstractProperty end
abstract type AbstractStructure end

const Molecule = Metagraph

function id(mol::Molecule)
    props(mol)[:id]
end

function position(mol::Molecule, atomid)
    props(mol, atomid)[:position]
end

function element(mol::Molecule, atomid)
    props(mol, atomid)[:element]
end

function labels(mol::Molecule, atomid)
    props(mol, atomid)[:labels]
end

# MetaGraph使ってるので分子と原子のpropはgpropとvpropでよいのでは？
#mutable struct Trajectory{T <: AbstractProperty}
#    snapshot::Vector{System}
#    prop::T
#end

mutable struct Position{T <: AbstractFloat} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

struct Topology
    g::Vector{Vector{MetaGraph}}
    lifetime::NamedTuple
end

mutable struct TimeSeries{T <: AbstractFloat}
    step::Vector{T}
    time::Vector{T}
    prop::Dict{HeierchyLabel, Any}
end

mutable struct EventSeries

end
"""

"""
mutable struct Trajectory{T <: AbstractFloat}
    # グラフと共に不変な情報はMetagraph内のハッシュテーブルに格納
    topology::EvecnSeries{Topology} #with atomid, charge, element, node label, edge label
    timeseries::TimeSeries{T}
    prop::Dict{Symbol, Any}
end

mutable struct System
    molecule::Vector{Molecule}
    box::BoundingBox
end



#@select {atoms | step ∈ [1:100], molid ∈ [1:10], atomid ∈ [1:100], label == "main chain"}
macro select_str(ex)
    str
    :(
        for i in 1:length($ex)
            println(string($ex[i]))
        end
    )
end

# 任意の時間，分子，ラベル，id，typeの原子に対する参照を取得できるようにする
# グラフ情報を保つ or not
# 任意の時間，分子，ラベル，id，typeの原子に対する参照を取得できるようにする
function select(traj, sym; step, mol, atom, label_condition)
    filter(x, view(traj, step))
end

function molecule_graph(traj; step, molid)
    
end


end