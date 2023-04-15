module HMD

using Graphs
using HDF5
using LinearAlgebra
using PeriodicTable
using Reexport
using SimpleWeightedGraphs
using StaticArrays
#using Unitful

import Base: getindex, firstindex, lastindex, setproperty!, iterate, length, precision, close, string, show, showerror
import Base: >, <, >=, <=, +, -, *, /, ==, position
import Base: position, time, contains, promote_type, promote_rule, similar
import Graphs: neighbors

# utility
export oblique_coord

# LabelHierarchy
export HLabel, LabelHierarchy, id, type

export
    # System constants
    Entire_System,

    # system core interface
    AbstractBbox,
    AbstractSystem,
    AbstractSystemType,
    AbstractTrajectory,
    dimension,
    precision,
    system_type,
    similar,
    show,
    natom,
    nbond,
    time,
    set_time!,
    topology,
    box,
    set_box!,
    all_elements,
    element,
    element,
    set_element!,
    set_elements!,
    all_positions,
    position,
    set_position!,
    set_position!,
    all_travels,
    travel,
    set_travel!,
    wrapped,
    wrap!,
    unwrap!,
    label2atom,

    # system label manipulation
    hierarchy_names,
    hierarchy,
    add_hierarchy!,
    remove_hierarchy!,
    all_labels,
    add_label!,
    add_labels!,
    count_label,
    add_relation!,
    add_relations!,
    insert_relations!,
    remove_label!,
    remove_relation!,
    contains,
    issuper,
    issub,
    super,
    sub,

    # system io interface
    AbstractFileFormat,
    close,

    # trajectory interface
    get_system,
    all_timesteps,
    get_timestep,
    length,
    add!,
    import_dynamic!,
    import_static!,
    latest_reaction,
    similar,
    similar_system,
    dimension,
    precision,
    system_type,
    wrapped,
    wrap!,
    unwrap!,
    add_snapshot!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    wrapped,

    # trajectory io interface
    add_snapshot!,
    import_dynamic!,
    import_static!,
    latest_reaction_step,
    get_timesteps,
    get_reactions,
    get_metadata,
    is_reaction,
    length,
    wrapped,

    # System highlevel interfaces
    is_atom,
    atom_label,
    add_atom!,
    add_atoms!,
    add_bond!,
    add_bonds!,
    bond_order,
    valence,
    neighbors,
    l2a,
    super_labels,
    sub_labels,
    hmdsave,
    read_system,

    # Trajectory highlevel interfaces
    length,
    getindex,
    firstindex,
    lastindex,
    iterate,
    similar,
    iterate,
    hmdsave,
    read_traj,
    snapshot

include("util.jl")

include("HierarchyLabels/HierarchyLabels.jl")
using .HierarchyLabels

include("interface/interface.jl")

include("DataTypes/DataTypes.jl")
@reexport using .DataTypes

#import .DataTypes: hmdsave, hmdread!
include("interface/system_highlevel.jl")
include("interface/trajectory_highlevel.jl")

#include("GenericIO/GenericIO.jl")
#using .GenericIO

end
