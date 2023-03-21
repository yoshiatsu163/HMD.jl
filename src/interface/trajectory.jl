function Base.getindex(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if !(0 < index <= length(traj))
        throw(BoundsError(traj, index))
    end

    replica = System{D, F, SysType}()
    # set properties that changes only at reaction
    rp = get_system(traj, latest_reaction(traj, index))
    replica.element = all_elements(rp) |> deepcopy
    replica.topology = topology(rp) |> deepcopy
    replica.hierarchy = deepcopy(rp.hierarchy)

    # set properties that changes at every step
    set_time!(replica, time(s))
    set_box!(replica, deepcopy(box(rp)))
    replica.position = all_positions(rp) |> deepcopy
    replica.travel = deepcopy(rp.travel)

    # others
    replica.wrapped = rp.wrapped

    return replica
end

#
function update_reader!(reader::System{D, F, Immutable}, traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    current = _get_system(traj, index)
    latest = latest_reaction(traj, index)[2]
    set_time!(reader, time(current))
    set_box!(reader, box(current))
    reader.position = all_positions(current)
    reader.element = all_elements(current)
    reader.travel = current.travel
    reader.wrapped = latest.wrapped
    reader.hierarchy = latest.hierarchy
    reader.topology = topology(latest)
    # TODO: props

    return nothing
end

#function Base.setproperty!(s::System{D, F, Immutable}, fieldname::Symbol) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
#    error("""This type $(typeof(s)) is intended to be read-only. If you want to mutate some data in trajectory, "s = traj[i]" makes a deepcopy. """)
#end

function Base.iterate(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    # ここでreaction pointをすべて調べ田植えでstateとして後続に渡せば毎回線型探索しなくても良い
    reaction = latest_reaction(traj, 1)
    index = 1
    reader = System{D, F, Immutable}()
    for field in fieldnames(typeof(s))
        setfield!(reader, field, getfield(s, field))
    end

    return (step=get_timestep(traj, index), reader=reader), index+1
    #return reader, timestep
end

function Base.iterate(traj::Trajectory{D, F, SysType}, state::Int64) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    index = state
    if index < length(traj)
        reader = System{D, F, Immutable}()
        return (step=get_timestep(traj, index), reader=reader), index+1
    else
        return nothing
    end
end

function wrap!(traj::AbstractTrajectory)
    for s in traj.systems
        wrap!(s)
    end

    return nothing
end

function unwrap!(traj::AbstractTrajectory)
    for s in traj.systems
        unwrap!(s)
    end

    return nothing
end

function wrapped(traj::AbstractTrajectory)
    is_wrapped = wrapped(traj.systems[1])
    @assert all(s -> wrapped(s)==is_wrapped, traj.systems)

    return is_wrapped
end


#####
##### Trajectory HDF5 interface
#####

function hmdsave(name::AbstractString, traj::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    file_handler = h5traj(name, "w")
    for reader in traj
        add_snapshot!(file_handler, reader.reader, step=reader.step)
    end
    close(file)
end

function read_traj(name::AbstractString)
    traj_file = h5traj(name, "r")

    D, F, SysType = get_metadata(traj_file)
    traj = Trajectory{D, F, SysType}()

    for reader in traj_file
        add!(traj, reader.reader, reader.step, change=DataType._is_reaction(reader.reader))
    end

    return traj
end

function Base.getindex(traj_file::H5traj, index::Integer)
    D, F, SysType = get_metadata(traj_file)
    s = System{D, F, SysType}()
    import_dynamic!(s, traj_file, index)
    latest_reaction = latest_reaction(traj_file, index)
    import_static!(s, traj_file, latest_reaction)

    return deepcopy(s)
end

function Base.iterate(traj_file::H5traj)
    index = 1

    reactions = get_reactions(traj_file)
    reaction_points = findall(==(true), reactions)
    timesteps = get_timesteps(traj_file)

    reader = System{D, F, Immutable}()
    import_static!(reader, traj_file, index)
    static_buffer = deepcopy(reader)
    import_dynamic!(reader, traj_file, index)

    return (step=step, reader=reader), (index+1, reaction_points, timesteps, static_buffer)
end

function Base.iterate(traj_file::H5traj, state::Tuple{Int64, Vector{Int64}, Vector{Int64}, System{D, F, SysType}}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    index, reaction_points, timesteps, static_buffer = state
    if !(0 < index < length(traj_file))
        throw(BoundsError(traj_file, index))
    elseif index == length(traj_file)
        return nothing
    end

    reader = System{D, F, Immutable}()
    if index ∈ reaction_points
        import_static!(reader, traj_file, index)
        static_buffer = deepcopy(reader)
    else
        reader.wrap = static_buffer.wrap
        reader.element = static_buffer.element
        reader.topology = static_buffer.topology
        reader.hierarchy = static_buffer.hierarchy
    end
    import_dynamic!(reader, traj_file, index)

    return (step=step, reader=reader), (index+1, reaction_points, timesteps, static_buffer)
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
