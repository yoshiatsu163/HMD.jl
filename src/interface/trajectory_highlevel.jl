function getindex(traj::AbstractTrajectory{D, F}, index::Integer) where {D, F<:AbstractFloat}
    if !(0 < index <= length(traj))
        throw(BoundsError(traj, index))
    end

    replica = System{D, F, Immutable}()
    # set properties that changes only at reaction
    rp = get_system(traj, latest_reaction(traj, index))
    replica.element = all_elements(rp) |> deepcopy
    replica.topology = topology(rp) |> deepcopy
    replica.hierarchy = deepcopy(rp.hierarchy)

    # set properties that changes at every step
    current = get_system(traj, index)
    set_time!(replica, time(replica))
    set_box!(replica, deepcopy(box(current)))
    replica.position = all_positions(current) |> deepcopy
    replica.travel = deepcopy(current.travel)

    # others
    replica.wrapped = current.wrapped

    return replica
end

function firstindex(traj::AbstractTrajectory)
    return 1
end

function lastindex(traj::AbstractTrajectory)
    return length(traj)
end

#function Base.setproperty!(s::System{D, F, Immutable}, fieldname::Symbol) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#    error("""This type $(typeof(s)) is intended to be read-only. If you want to mutate some data in trajectory, "s = traj[i]" makes a deepcopy. """)
#end

function Base.iterate(traj::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    index = 1
    reader = System{D, F, Immutable}()
    DataTypes.import_dynamic!(reader, traj, index)
    DataTypes.import_static!(reader, traj, index)

    return (step=get_timestep(traj, index), reader=reader), index+1
end

function Base.iterate(traj::AbstractTrajectory{D, F}, state::Int64) where {D, F<:AbstractFloat}
    index = state
    if index <= length(traj)
        reader = System{D, F, Immutable}()
        rp = latest_reaction(traj, index)
        DataTypes.import_static!(reader, traj, rp)
        DataTypes.import_dynamic!(reader, traj, index)
        return (step=get_timestep(traj, index), reader=reader), index+1
    else
        return nothing
    end
end

#####
##### Trajectory HDF5 interface
#####

function hmdsave(name::AbstractString, traj::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    file_handler = h5traj(name, "w")
    index = 1
    for reader in traj
        add_snapshot!(file_handler, reader.reader, reader.step; reaction=is_reaction(traj, index), unsafe=true)
        index += 1
    end
    close(file_handler)
end

function read_traj(name::AbstractString, template::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    traj_file = h5traj(name, "r")

    D_file, F_file, SysType_file = get_metadata(traj_file)
    if (D_file, F_file) != (D, F)
        error("Trajectory file type ($D_file, $F_file) is not compatible with the template type ($D, $F).")
    end

    traj = similar(template)
    timesteps = get_timesteps(traj_file)
    reaction_points = get_reactions(traj_file)
    for (step, index) in zip(timesteps, 1:length(traj_file))
        s = snapshot(traj_file, index, similar_system(template))
        add!(traj, s, step; reaction=(step ∈ reaction_points))
    end

    return traj
end

function snapshot(traj_file::H5traj, index::Integer, template::AbstractSystem{D, F}) where {D, F<:AbstractFloat}
    D_file, F_file, SysType_file = get_metadata(traj_file)
    if (D_file, F_file) != (D, F)
        error("Trajectory file $(name) is not compatible with the template $(template).")
    end
    step = get_timesteps(traj_file)[index]

    s = similar(template)
    DataTypes.import_dynamic!(s, traj_file; step=step)
    latest_react = latest_reaction_step(traj_file, step)
    DataTypes.import_static!(s, traj_file, step=latest_react)

    return s
end

function Base.iterate(traj_file::H5traj)
    index = 1

    # timestep at reaction (not index!)
    reaction_steps = get_reactions(traj_file)
    timesteps = get_timesteps(traj_file)

    D, F, stub = get_metadata(traj_file)
    reader = System{D, F, Immutable}()
    DataTypes.import_static!(reader, traj_file, step=timesteps[index])
    static_cache = deepcopy(reader)
    DataTypes.import_dynamic!(reader, traj_file, step=timesteps[index])

    return (step=timesteps[index], reader=reader), (index+1, reaction_steps, timesteps, static_cache)
end

function Base.iterate(traj_file::H5traj, state::Tuple{Int64, Vector{Int64}, Vector{Int64}, System{D, F, Immutable}}) where {D, F<:AbstractFloat}
    index, reaction_steps, timesteps, static_cache = state
    if index > length(traj_file)
        return nothing
    end

    reader = System{D, F, Immutable}()
    if timesteps[index] ∈ reaction_steps
        DataTypes.import_static!(reader, traj_file; step=timesteps[index])
        static_cache = deepcopy(reader)
    else
        reader.wrapped = static_cache.wrapped
        reader.element = static_cache.element
        reader.topology = static_cache.topology
        reader.hierarchy = static_cache.hierarchy
    end
    DataTypes.import_dynamic!(reader, traj_file; step=timesteps[index])

    return (step=timesteps[index], reader=reader), (index+1, reaction_steps, timesteps, static_cache)
end

# getindex?
#function slice(traj::AbstractTrajectory, index::Integer)
#
#end
#
#function slice(traj::AbstractTrajectory, time::AbstractFloat)
#
#end
#
#function to_system(traj::AbstractTrajectory)
#
#end
#
#function nearest_slice(traj::AbstractTrajectory, time::AbstractFloat)
#
#end
#
#function Base.push!(traj::AbstractTrajectory{D, F, SysType}, s::AbstractSystem{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
#
#function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
#
#function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#
#end
