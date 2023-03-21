export AbstractTrajectory, Immutable, Trajectory
export all_timesteps, get_timestep, is_reaction, get_system
export latest_reaction, add!, update_reader!, add!
export setproperty!, iterate, getindex, length
#export hmdsave, hmdread, add_snapshot!

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

function is_reaction(traj::Trajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return traj.is_reaction
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

    return findlast(==(true), view(is_reaction(traj), 1:index))
end


#####
##### Trajectory HDF5 types
#####

function add_snapshot!(file_handler::H5traj, reader::System{D, F, SysType}, step) where{D, F<:AbstractFloat, SysType<:AbstractSystemType}
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
    else
        _error_chk(file, reader; mode=mode)
    end

    # trajectory specific properties
    file["reactions/$step"] = _is_reaction(reader)
    snap = H5system(file["snapshots/$step"])
    hmdsave(snap, reader)
    
    return nothing
end

function snapshot(file_handler::H5traj, index::Integer)
    file = get_file(file_handler)
    timesteps = keys(file["snapshots"])

    step = parse(Int64, timesteps[index])
    snapshot = h5system(file["snapshots/$step"]) |> read_system
    if _is_reaction(snapshot)
        return read_system(snap)
    else
        reactions = keys(file["reactions"])
        i = findlast(<(step), reactions)
        latest_reaction = 
    end
end

function nstep(file_handler::H5traj)
    file = get_file(file_handler)
    keys(file) |> filter(x->occursin(r"\d+", x)) |> length
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

function _is_reaction(s::AbstractSystem)
    return all_positions(s) |> isempty
end