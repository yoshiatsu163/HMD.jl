module HMD

using Graphs
using JLD2
using LinearAlgebra
using MetaGraphs
using PeriodicTable
using StaticArrays
using Unitful

export Position, BoundingBox, AbstractSystem, System, Label
export time, set_time!, natom, topology, box, set_box!
export all_element, element, set_element!
export all_positions, position, set_position!
export hierarchy_names, hierarchy, add_hierarchy!, remove_hierarchy!, merge_hierarchy!
export prop_names, props, prop, labels_in_prop, add_prop!, set_prop!
export labels, add_label!, add_relation!, insert_relation!, remove_label!, remove_relation!
export Id, Category
export contains, has_relation, issuper, issub, super, sub
export hmdsave, hmdread

include("DataTypes/DataTypes.jl")
using .DataTypes
include("interface.jl")

include("GenericIO/GenericIO.jl")
using .GenericIO

end
