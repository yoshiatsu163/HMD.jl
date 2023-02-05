module HMD

using Graphs
using JLD2
using LinearAlgebra
using MetaGraphs
using PeriodicTable
using SimpleWeightedGraphs
using StaticArrays
using Unitful

import Base: getindex
import Base: >, <, >=, <=, +, -, *, /, ==, string, show
import Base: position, time, contains, show
import MetaGraphs: set_prop!, props
import Graphs: neighbors

export Position, BoundingBox, AbstractSystem, System, HLabel, LabelHierarchy
export time, set_time!, natom, topology, box, set_box!, add_atom!, add_bond!, atom_label, l2a, is_atom
export all_elements, element, set_element!
export all_positions, position, set_position!
export hierarchy_names, hierarchy, add_hierarchy!, remove_hierarchy!, merge_hierarchy!
export prop_names, props, prop, labels_in_prop, add_prop!, set_prop!
export labels, add_label!, count_label, add_relation!, insert_relation!, remove_label!, remove_relation!
export Id, Category, Entire_System
export id, type, ==
export contains, has_relation, issuper, issub, super, sub
export hmdsave, hmdread
export valence, bond_order, neighbors

include("DataTypes/DataTypes.jl")
using .DataTypes
include("interface.jl")

include("GenericIO/GenericIO.jl")
using .GenericIO

end
