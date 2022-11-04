module HMDs

using Graphs
using MetaGraphs
using PeriodicTable
using Unitful
using StaticArrays
using Pipe

include("GenericIO/GenericIO.jl")
using .GenericIO

include("HeierchyLabels.jl")
using .HeierchyLabels

# 人間絶叫多数、騒乱
mutable struct HMD{T <: Integer, U <: AbstractFloat}
    graph::MetaGraph{T, U}
    graph
end
"""
グラフ間の包含関係や頻度情報はGraphManipulationで扱う
    グラフが重複することもあり得る
    グラフ間の包含関係はインスタンスが無ければわからない
        従ってグラフラベルの段階では包含関係を定義しない
ここでは単にグラフと原子の所属関係を保存する
    ・グラフ名->ノードid[] in gprop
    ・グラフ名 in each vprop
Labelの設計
グラフのイテレータも必要
    余) イテレータの定義はオートマトンの定義そのもの？
グラフラベルへのアクセスパターン
    some subgraph -> {原子, property, subgraph, supergraph}
"""

#include("TrajectoryContainer/")




end
