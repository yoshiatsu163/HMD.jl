module HMD

using Graphs
using JLD2
using LinearAlgebra
#using MetaGraphs
using PeriodicTable
using SimpleWeightedGraphs
using StaticArrays
using Unitful

import Base: getindex, setproperty!, iterate, length, precision
import Base: >, <, >=, <=, +, -, *, /, ==, string, show
import Base: position, time, contains, show, promote_type, promote_rule
#import MetaGraphs: set_prop!, props
import Graphs: neighbors

export id, type, ==, promote_rule, promote_type, position, time, contains, show, name, StaticString
export >, <, >=, <=, +, -, *, /, ==, string, show, convert, getindex, convert, setproperty!, iterate

export AbstractSystemType, GeneralSystem
export Position, BoundingBox, AbstractSystem, System, HLabel, H_Label, LabelHierarchy
export time, set_time!, natom, nbond, topology, box, set_box!, dimension, precision, add_atom!, add_atoms!, add_bond!, add_bonds!, atom_label, l2a, is_atom
export all_elements, element, _add_element!, set_element!
export all_positions, position, _add_position!, set_position!, all_travels, travel, set_travel!, wrapped, wrap!, unwrap!
export hierarchy_names, hierarchy, add_hierarchy!, remove_hierarchy!, merge_hierarchy!
export prop_names, props, prop, labels_in_prop, add_prop!, set_prop!
export labels, add_label!, add_labels!, count_label, add_relation!, add_relations!, insert_relation!, remove_label!, remove_relation!
export Id, Category, Entire_System
export id, type, ==
export contains, has_relation, issuper, issub, super, sub
export hmdsave, hmdread
export dimension, valence, bond_order, neighbors, all_labels, super_labels, sub_labels, wrap, atom_mass

export AbstractTrajectory, Immutable, Trajectory
export all_times, all_timesteps, get_timestep, timestep2time, timestep2index, change_points
export latest_changepoint, add!, update_reader!
export setproperty!, iterate, getindex, length

include("DataTypes/DataTypes.jl")
using .DataTypes
include("interface.jl")

#include("GenericIO/GenericIO.jl")
#using .GenericIO

end
