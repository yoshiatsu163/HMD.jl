abstract type AbstractSystem{D} end



# data APIs ================================================================================================

"""
    subspace(condition)

A interace for temporal, spacial, and structural slicing of trajectory.

⋅⋅⋅
# Arguments
- 'condition': condition to specify region in a trajectory.
⋅⋅⋅

"""
function subspace end


"""
    read_prop(condition)

A interace for reading properties of subspace.

⋅⋅⋅
# Arguments
- 'condition': condition to specify region in a trajectory.
⋅⋅⋅

"""
function read_prop end


"""
write(data)

A interace for writing properties, BCs, metadata into subspace.

⋅⋅⋅
# Arguments
- 'data': Some data which can be wrote to a trajectory without changing the chemical structure in it.
⋅⋅⋅

"""
function write end


"""
⊕(rhs, lhs)

A interace for combining chemical stuctures.

⋅⋅⋅
# Arguments
- 'rhs, lhs': Right and left hand side subspace object whose time lengh == 1.
⋅⋅⋅

"""
function ⊕ end


"""
metadata(traj)

A interace for reading metadata.

⋅⋅⋅
# Arguments
- 'traj': A trajectory object.
⋅⋅⋅


"""
function metadata end

"""
BC(traj)

A interace for reading boundary condition.

⋅⋅⋅
# Arguments
- 'traj': A trajectory object.
⋅⋅⋅

"""
function BC end


"""
new(args)

A interace for contstructing a trajectory.

⋅⋅⋅
# Arguments
- 'args': Arguments needed to contstruct a trajectory object.
⋅⋅⋅

"""
function new end


"""
heierchy(traj)

A interace for reading subgraph heierchy dictionary: name -> heierchy tree.

⋅⋅⋅
# Arguments
- 'traj': A trajectory object.
⋅⋅⋅

"""
function heierchy end


"""
bounding_box(traj)

A interace for reading subgraph heierchy dictionary: name -> heierchy tree.

⋅⋅⋅
# Arguments
- 'traj': A trajectory object.
⋅⋅⋅

"""
function bounding_box(sys <: AbstractSystem) end


# disc IO APIs ================================================================================================
#"""
#heierchy(traj)
#
#A interace for combining chemical stuctures.
#
#⋅⋅⋅
## Arguments
#- 'traj': A trajectory object.
#⋅⋅⋅
#
#"""
#function  end
#
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end
#"""
#"""
#function  end

