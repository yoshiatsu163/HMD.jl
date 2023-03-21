export AbstractTrajectory, Immutable, Trajectory
export all_timesteps, get_timestep, is_reaction, get_system
export latest_reaction, add!, update_reader!, add!
export setproperty!, iterate, getindex, length
export add_snapshot!, get_timesteps, get_reactions, get_metadata

abstract type AbstractTrajectory{D, F<:AbstractFloat} end
struct Immutable <: AbstractSystemType end

#use case
# 原子位置情報の追跡
# 時系列propの追跡

#設計
# 1. buffer::systemを与えて時系列データのポインタだけ差し替え
# 2. Systemのbase関数を多重ディスパッチでtrajctoryに拡張 read_onlyならいける?
# -> 2を採用(2の内部で結局1と同じことを実行する)


# 変化の時間スケールに差がある変数を分ける
Base.@kwdef mutable struct Trajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractTrajectory{D, F}
    systems::Vector{System{D, F, SysType}} = System{D, F, SysType}[]
    is_reaction::Vector{Bool} = Vector{Bool}(undef, 0)
    timesteps::Vector{Int64} = Vector{Int64}(undef, 0)
end

function Trajectory(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return Trajectory{D, F, SysType}([s], [true], [1])
end

function is_reaction(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.is_reaction[index]
end

function get_system(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.systems[index]
end

function all_timesteps(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.timesteps
end

function get_timestep(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return all_timesteps(traj)[index]
end

function Base.length(traj::Trajectory{D, F, SysType}{D, F}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return length(traj.systems)
end

function add!(traj::Trajectory{D, F, SysType}, s::System{D, F, SysType}, timestep::Integer; change=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    @assert length(traj.systems) == length(traj.timesteps)
    @assert length(traj.is_reaction) == length(traj)

    push!(traj.timesteps, timestep)
    push!(traj.is_reaction, change)
    if change
        push!(traj.systems, deepcopy(s))
    else
        replica = System{D, F, SysType}()
        set_box!(replica, deepcopy(box(s)))
        set_time!(replica, time(s))
        replica.position = all_positions(s) |> deepcopy
        replica.travel = deepcopy(s.travel)
        replica.wrap = s.wrap
        #replica.props = deepcopy(s.props)
        push!(traj.systems, replica)
    end

    return nothing
end

function latest_reaction(traj::Trajectory{D, F, SysType}, index::Integer) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if index < 1
        error("index must be positive. ")
    end

    return findlast(==(true), view(traj.is_reaction, 1:index))
end

function is_reaction(s::AbstractSystem)
    return all_positions(s) |> isempty
end

#####
##### Trajectory HDF5 types
#####

function add_snapshot!(file_handler::H5traj, reader::System{D, F, SysType}, step; unsafe=false) where{D, F<:AbstractFloat, SysType<:AbstractSystemType}
    file = get_file(file_handler)
    # construction
    if keys(file) |> isempty
        file["infotrype"] = "Trajectory"
        file["dimension"] = dimension(reader)
        file["precision"] = precision(reader) |> string
        file["system_type"] = system_type(reader) |> string
        file["wrapped"] = wrapped(reader)
        create_group(file, "snapshots")
        create_group(file, "reactions")
    elseif !unsafe
        _error_chk(file, reader; mode=mode)
    end

    # trajectory specific properties
    file["reactions/$step"] = is_reaction(reader)
    snap = H5system(file["snapshots/$step"])
    hmdsave(snap, reader)

    return nothing
end

function import_dynamic!(reader::System{D, F, Immutable}, traj_file::H5traj, index::Integer; timesteps::Vector{Int64}=Int64[]) where {D, F<:AbstractFloat}
    if isempty(timesteps)
        timesteps = get_timesteps(traj_file)
    end

    file = get_file(traj_file)
    snap = h5system(file["snapshots/$(timesteps[index])"], "r")
    D_file, F_file, SysType_file = get_metadata(snap)
    if (D, F) != (D_file, F_file)
        error("Dimension and precision of the reader and the trajectory file are different. ")
    end
    import_dynamic!(reader, snap)

    return nothing
end

function import_static!(reader::System{D, F, Immutable}, traj_file::H5traj, index::Integer; timesteps::Vector{Int64}=Int64[]) where {D, F<:AbstractFloat}
    if isempty(timesteps)
        timesteps = get_timesteps(traj_file)
    end

    file = get_file(traj_file)
    snap = h5system(file["snapshots/$(timesteps[index])"], "r")
    D_file, F_file, SysType_file = get_metadata(snap)
    if (D, F) != (D_file, F_file)
        error("Dimension and precision of the reader and the trajectory file are different. ")
    end
    import_static!(reader, snap)

    return nothing
end

function latest_reaction(traj_file::H5traj, index::Integer)
    reactions = get_reactions(traj_file)
    return findlast(==(true), view(reactions, 1:index))
end

function get_timesteps(traj_file::H5traj)
    file = get_file(traj_file)
    return keys(file["snapshots"])
end

function get_reactions(traj_file::H5traj)
    file = get_file(traj_file)
    timestep = keys(file["reactions"])
    reactions = Vector{Bool}(undef, length(timestep))
    for step in timestep
        reactions[step] = read(file["reactions/$step"])
    end

    return reactions
end

function get_metadata(traj_file::H5traj)
    file = get_file(traj_file)
    D = read(file, "dimension")
    F = read(file, "precision") |> Symbol |> eval
    SysType = read(file, "system_type") |> Symbol |> eval

    return D, F, SysType
end

function _error_chk(file, reader::AbstractSystem{D, F}; mode) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if read(file, "infotrype") != "Trajectory"
        error("The file is not a trajectory file. ")
    elseif read(file, "dimension") != dimension(reader)
        error("The dimension of the system is different from the dimension of the trajectory. ")
    elseif read(file, "precision") != string(precision(reader))
        error("The precision of the system is different from the precision of the trajectory. ")
    elseif read(file, "system_type") != string(system_type(reader))
        error("The system type of the system is different from the system type of the trajectory. ")
    elseif read(file, "mode") != mode
        error("mode mistmatch")
    end

    if wrapped(reader)
        error("currently supports only unwrapped coordinates. ")
    end
end

function is_reaction(traj_file::H5traj, index::Integer)
    file = get_file(traj_file)
    return read(file, "reactions/$index")
end

function Base.length(traj_file::H5traj)
    file = get_file(traj_file)
    return length(keys(file["snapshots"]))
end
