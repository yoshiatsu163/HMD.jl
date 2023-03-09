#####
##### Id definition
#####

struct Id{T}
    id::Int64
end

function Int(i::Id)
    i.id
end

function Base.convert(::Type{<:Integer}, i::Id)
    Int64(i.id)
end

Base.promote_rule(::Type{<:Integer}, ::Type{Id{<:Any}}) = Id
#Base.promote_rule(::Type{Id}, ::Type{<:Integer}) = Id

for op in [:(+), :(-), :(*), :(/)]
    @eval function $op(lhs::Id{T}, rhs::Id{T}) where {T}
        Id{T}($op(lhs.id, rhs.id))
    end
    @eval function $op(lhs::Id{T}, rhs::Int64) where {T}
        Id{T}($op(lhs.id, rhs))
    end
    @eval function $op(lhs::Int64, rhs::Id{T}) where {T}
        Id{T}($op(lhs, rhs.id))
    end
end

for op in [:(>), :(<), :(>=), :(<=), :(==)]
    @eval function $op(lhs::Id{T}, rhs::Id{T}) where {T}
        $op(lhs.id, rhs.id)
    end
    @eval function $op(lhs::Id{T}, rhs::Int64) where {T}
        $op(lhs, rhs.id)
    end
    @eval function $op(lhs::Int64, rhs::Id{T}) where {T}
        $op(lhs.id, rhs)
    end
end

#####
##### Category definition
#####

struct Category{T}# <: AbstractString
    str::String
end

function string(category::Category)
    category.str
end

function ==(lhs::Category, rhs::Category)
    string(lhs) == string(rhs)
end

function Base.length(category::Category)
    return length(string(category))
end
