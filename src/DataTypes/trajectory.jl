abstract type AbstractTrajectory end

#use case
# 原子位置情報の追跡
# 時系列propの追跡

#設計
# 1. buffer::systemを与えて時系列データのポインタだけ差し替え
# 2. Systemのbase関数を多重ディスパッチでtrajctoryに拡張 read_onlyならいける?
# -> 2を採用(2の内部で結局1と同じことを実行する)
print_to_string

export contains
export all_elements, element
export time, natom, nbond, topology, box, dimension, show
export all_positions, position, travel, wrapped
export hierarchy_names, hierarchy
export all_labels, count_label, has_relation, issuper, issub, super, sub
export prop_names, props, prop, labels_in_prop


# 変化の時間スケールに差がある変数を分ける
mutable struct Trajectory{S<:AbstractSystem{D, F}} <: AbstractTrajectory{D, F}
    systems::SortedDict{Int64, S}
    step2time::SortedDict{Int64, F}
    change_points::Vector{Int64}
end

function change_points(traj::AbstractTrajectory)
    return traj.change_points
end

function _get_system(traj::AbstractTrajectory, timestep::Integer)
    return traj.systems[timestep]
end

function time(traj::AbstractTrajectory)
    
end

function latest_changepoint(traj::AbstractTrajectory, current::Integer)
    index = findfirst(i -> i > current, change_points(traj))
    timestep = change_points(traj)[index]

    return _get_system(traj, timestep)
end

function Base.getindex(traj::AbstractTrajectory, timestep::Integer)
    s = latest_changepoint(traj, timestep) |> deepcopy
    current = _get_system(traj, timestep)
    set_time!(s, time(current))
    set_box!(s, deepcopy(box(current)))
    s.position = all_positions(current) |> deepcopy
    s.travel = travel(current) |> deepcopy
    s.props = current.props

    return s
end

function Base.iterate(traj::AbstractTrajectory)
    return 1, time(traj, 1), latest_changepoint(traj, 1)
end

function Base.iterate(traj::AbstractTrajectory, state::Integer)
    
end

# getindex?
function slice(traj::AbstractTrajectory, index::Integer)
    
end

function slice(traj::AbstractTrajectory, time::AbstractFloat)
    
end

function to_system(traj::AbstractTrajectory)
    
end

function nearest_slice(traj::AbstractTrajectory, time::AbstractFloat)
    
end

function Base.push!(traj::AbstractTrajectory{D, F, SysType}, s::AbstractSystem{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    
end

function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    
end

function Base.append!(addend::AbstractTrajectory{D, F, SysType}, augend::AbstractTrajectory{D, F, SysType}) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    
end
