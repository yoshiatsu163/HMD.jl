export AbstractTrajectory, Immutable, Trajectory
export all_times, all_timesteps, get_timestep, timestep2time, timestep2index, change_points
export latest_changepoint, add!, update_reader!
export setproperty!, iterate, getindex, length

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
mutable struct Trajectory{D, F<:AbstractFloat, SysType<:AbstractSystemType} <: AbstractTrajectory{D, F}
    systems::Vector{System{D, F, SysType}}
    step2time::Vector{F}
    change_points::Vector{Int64}
    timesteps::Vector{Int64}
end

function Trajectory(s::System{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return Trajectory{D, F, SysType}([s], [0.0], [1], [1])
end

function change_points(traj::AbstractTrajectory)
    return traj.change_points
end

function _get_system(traj::AbstractTrajectory, index::Integer)
    return traj.systems[index]
end

function all_times(traj::AbstractTrajectory)
    return traj.step2time
end

function all_timesteps(traj::AbstractTrajectory)
    return traj.timesteps
end

function get_timestep(traj::AbstractTrajectory, index::Integer)
    return all_timesteps(traj)[index]
end

function time(traj::AbstractTrajectory, index::Integer)
    return all_times(traj)[index]
end

function timestep2time(traj::AbstractTrajectory, timestep::Integer)
    index = timestep2index(traj, timestep)
    return time(traj, index)
end

function timestep2index(traj::AbstractTrajectory, timestep::Integer)
    return findfirst(==(timestep), all_timesteps(traj))
end

function Base.length(traj::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    return length(traj.systems)
end

function add!(traj::Trajectory{D, F, SysType}, s::System{D, F, SysType}, timestep::Integer; change=false) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    @assert length(traj.systems) == length(traj.step2time) == length(traj.timesteps)
    @assert length(traj.change_points) <= length(traj)

    push!(traj.step2time, time(s))
    push!(traj.timesteps, timestep)
    if change
        push!(traj.systems, deepcopy(s))
        push!(traj.change_points, length(traj.systems))
    else
        replica = deepcopy(s)
        empty!(replica.topology)
        empty!(replica.hierarchy)
        empty!(replica.elements)
        push!(traj.systems, replica)
    end
end

function latest_changepoint(traj::AbstractTrajectory, index::Integer)
    min_index, max_index = extrema(change_points(traj))
    if index < min_index
        error("given index $(index) is smaller than the first change point $min_index. ")
    end

    cp = change_points(traj)
    if length(cp) < 2
        return time(traj, cp[1]), _get_system(traj, cp[1])
    elseif max_index <= index
        return time(traj, cp[end]), _get_system(traj, cp[end])
    elseif index ∈ cp
        return time(traj, index), _get_system(traj, index)
    else
        latest_index = findfirst(index_at_change -> (index_at_change > index), cp) - 1
        return time(traj, cp[latest_index]), _get_system(traj, cp[latest_index])
    end
end

function Base.getindex(traj::AbstractTrajectory{D, F}, index::Integer) where {D, F<:AbstractFloat}
    if !(0 < index <= length(traj))
        throw(BoundsError(traj, index))
    end

    s = latest_changepoint(traj, index)[2]
    replica = deepcopy(s)

    if index ∈ change_points(traj)
        return replica
    else
        current = _get_system(traj, index)
        set_time!(replica, time(current))
        set_box!(replica, deepcopy(box(current)))
        replica.position = all_positions(current) |> deepcopy
        replica.travel = deepcopy(current.travel)
        replica.props = deepcopy(current.props)
        return replica
    end
end

function update_reader!(reader::System{D, F, Immutable}, traj::AbstractTrajectory{D, F}, index::Integer) where {D, F<:AbstractFloat}
    current = _get_system(traj, index)
    latest = latest_changepoint(traj, index)[2]
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

#function Base.setproperty!(s::System{D, F, Immutable}, fieldname::Symbol) where {D, F<:AbstractFloat}
#    error("""This type $(typeof(s)) is intended to be read-only. If you want to mutate some data in trajectory, "s = traj[i]" makes a deepcopy. """)
#end

function Base.iterate(traj::AbstractTrajectory{D, F}) where {D, F<:AbstractFloat}
    s = latest_changepoint(traj, 1)[2]
    index = 1
    reader = System{D, F, Immutable}()
    for field in fieldnames(typeof(s))
        setfield!(reader, field, getfield(s, field))
    end

    return (step=get_timestep(traj, index), time=time(traj, index), reader=reader), (reader, index+1)
    #return reader, timestep
end

# 並列化時に競合の可能性
function Base.iterate(traj::AbstractTrajectory{D, F}, state::Tuple{System{D, F, Immutable}, Int64}) where {D, F<:AbstractFloat}
    reader, index = state
    if index < length(traj)
        update_reader!(reader, traj, index)
        return (step=get_timestep(traj, index), time=time(traj, index), reader=reader), (reader, index+1)
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
