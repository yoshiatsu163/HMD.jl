"""
    NotImplementedError{M}(m)
`Exception` thrown when a method from the `AbstractSystem` or `AbstractTrajectory` interface
is not implemented by a given AbstractSystem or Trajectory type.
"""
struct NotImplementedError{M} <: Exception
    m::M
    NotImplementedError(m::M) where {M} = new{M}(m)
end

function Base.showerror(io::IO, ie::NotImplementedError)
    return print(io, "method $(ie.m) not implemented.")
end

_NI(m) = throw(NotImplementedError(m))


abstract type AbstractSystem{D, F<:AbstractFloat} end
abstract type AbstractSystemType end
abstract type AbstractTrajectory{D, F<:AbstractFloat} end
abstract type AbstractBbox{D, F<:AbstractFloat} end

# system core interface
dimension(s::AbstractSystem) = _NI("dimension")
precision(s::AbstractSystem) = _NI("precision")
system_type(s::AbstractSystem) = _NI("system_type")
similar(s::AbstractSystem) = _NI("similar")
show(io::IO, ::MIME"text/plain", s::AbstractSystem) = _NI("show")
natom(s::AbstractSystem) = _NI("natom")
nbond(s::AbstractSystem) = _NI("nbond")
time(s::AbstractSystem) = _NI("time")
set_time!(s::AbstractSystem, time::AbstractFloat) = _NI("set_time!")
topology(s::AbstractSystem) = _NI("topology")
box(s::AbstractSystem) = _NI("box")
set_box!(s::AbstractSystem, box::AbstractBbox) = _NI("set_box!")
all_elements(s::AbstractSystem) = _NI("all_elements")
element(s::AbstractSystem, i::Integer) = _NI("element")
element(s::AbstractSystem, label::HLabel) = _NI("element")
set_element!(s::AbstractSystem, atom_id::Integer, ename::AbstractString) = _NI("set_element!")
set_elements!(s::AbstractSystem, atom_ids::AbstractVector{<:Integer}, enames::AbstractVector{<:AbstractString}) = _NI("set_elements!")
all_positions(s::AbstractSystem) = _NI("all_positions")
position(s::AbstractSystem, atom_id::Integer) = _NI("position")
set_position!(s::AbstractSystem, atom_id::Integer, x::AbstractVector{<:AbstractFloat}) = _NI("set_position!")
set_position!(s::AbstractSystem, label::HLabel, x::AbstractVector{<:AbstractFloat}) = _NI("set_position!")
all_travels(s::AbstractSystem) = _NI("all_travels")
travel(s::AbstractSystem, atom_id::Integer) = _NI("travel")
set_travel!(s::AbstractSystem, atom_id::Integer, n::AbstractVector{<:Integer}) = _NI("set_travel!")
wrapped(s::AbstractSystem) = _NI("wrapped")
wrap!(s::AbstractSystem) = _NI("wrap!")
unwrap!(s::AbstractSystem) = _NI("unwrap!")

# system label manipulation
hierarchy_names(s::AbstractSystem) = _NI("hierarchy_names")
hierarchy(s::AbstractSystem, name::AbstractString) = _NI("hierarchy")
add_hierarchy!(s::AbstractSystem, name::AbstractString) = _NI("add_hierarchy!")
remove_hierarchy!(s::AbstractSystem, name::AbstractString) = _NI("remove_hierarchy!")
all_labels(s::AbstractSystem, hname::AbstractString) = _NI("all_labels")
all_labels(s::AbstractSystem, hname::AbstractString, label_type::AbstractString) = _NI("all_labels")
add_label!(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("add_label!")
add_label!(s::AbstractSystem, hname::AbstractString, label_type::AbstractString) = _NI("add_label!")
add_labels!(s::AbstractSystem, hname::AbstractString, labels::AbstractVector{<:HLabel}) = _NI("add_labels!")
count_label(s::AbstractSystem, hname::AbstractString, label_type::String) = _NI("count_label")
add_relation!(s::AbstractSystem, hname::AbstractString; super::HLabel, sub::HLabel) = _NI("add_relation!")
add_relations!(s::AbstractSystem, hname::AbstractString; super::AbstractVector{<:HLabel}, sub::AbstractVector{<:HLabel}) = _NI("add_relations!")
insert_relations!(s::AbstractSystem, hname::AbstractString, index::Integer; super::AbstractVector{<:HLabel}, sub::AbstractVector{<:HLabel}) = _NI("insert_relations!")
remove_label!(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("remove_label!")
remove_relation!(s::AbstractSystem, hname::AbstractString, super::HLabel, sub::HLabel) = _NI("remove_relation!")
contains(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("contains")
issuper(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel) = _NI("issuper")
issub(s::AbstractSystem, hname::AbstractString, label1::HLabel, label2::HLabel) = _NI("issub")
super(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("super")
sub(s::AbstractSystem, hname::AbstractString, label::HLabel) = _NI("sub")

# property interface


# system_io interface
abstract type AbstractFileFormat end
close(file_handler::AbstractFileFormat) = _NI("close")
# TODO Abstract File Formatを元にsystemとtrajectoryの構成をつくる

# trajectory interface
is_reaction(s::AbstractTrajectory, index::Integer) = _NI("is_reaction")
is_reaction(s::AbstractSystem) = _NI("is_reaction")
get_system(s::AbstractTrajectory, index::Integer) = _NI("get_system")
all_timesteps(traj::AbstractTrajectory) = _NI("all_timesteps")
get_timestep(traj::AbstractTrajectory, index::Integer) = _NI("get_timestep")
length(traj::AbstractTrajectory) = _NI("length")
add!(traj::AbstractTrajectory, s::AbstractSystem, timestep::Integer; reaction=false) = _NI("add!")
latest_reaction(traj::AbstractTrajectory, index::Integer) = _NI("latest_reaction")
similar(traj::AbstractTrajectory) = _NI("similar")
similar_system(traj::AbstractTrajectory) = _NI("System")
dimension(traj::AbstractTrajectory) = _NI("dimension")
precision(traj::AbstractTrajectory) = _NI("precision")
system_type(traj::AbstractTrajectory) = _NI("system_type")
wrapped(traj::AbstractTrajectory) = _NI("wrapped")
wrap!(traj::AbstractTrajectory) = _NI("wrap!")
unwrap!(traj::AbstractTrajectory) = _NI("unwrap!")

# trajectory io interface
add_snapshot!(file_handler::AbstractFileFormat, s::AbstractSystem, step::Int64; reaction::Bool=false, unsafe::Bool=false) = _NI("add_snapshot!")
import_dynamic!(reader::AbstractSystem, traj_file::AbstractFileFormat; index::Int64=typemin(Int64), step::Int64=typemin(Int64)) = _NI("import_dynamic!")
latest_reaction_step(traj_file::AbstractFileFormat, current_step::Int64) = _NI("latest_reaction_step")
get_timesteps(traj_file::AbstractFileFormat) = _NI("get_timesteps")
get_reactions(traj_file::AbstractFileFormat) = _NI("get_reactions")
get_metadata(traj_file::AbstractFileFormat) = _NI("get_metadata")
is_reaction(traj_file::AbstractFileFormat, index::Int64) = _NI("is_reaction")
length(traj_file::AbstractFileFormat) = _NI("length")
wrapped(traj_file::AbstractFileFormat) = _NI("wrapped")
