#####
##### StaticString definition
#####

struct StaticString{N}
    str::NTuple{N, UInt8}
end

function StaticString(str::AbstractString)
    sstr = Tuple(UInt8(c) for c in str)
    return StaticString{length(sstr)}(sstr)
end

function string(sstr::StaticString)
    return Char.(sstr.str) |> join
end

function Base.show(io::IO, ::MIME"text/plain", sstr::StaticString)
    print(io, "$(string(sstr))")
end

function Base.print_to_string(sstr::StaticString)
    return "$(string(sstr))"
end

function *(lhs::StaticString, rhs::StaticString)
    return StaticString(string(lhs) * string(rhs))
end

function Base.iterate(sstr::StaticString)
    return if isempty(sstr.str)
        nothing
    else
        sstr.str[1], 2
    end
end

function Base.iterate(sstr::StaticString, state::Integer)
    if state < length(sstr.str)
        return sstr.str[state], state + 1
    else
        return nothing
    end
end

function Base.length(sstr::StaticString)
    return length(sstr.str)
end

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

struct Category{T}
    str::StaticString
end

function Category{T}(str::AbstractString) where {T}
    return Category{T}(StaticString(str))
end

function Category{T}(vec::NTuple{N, UInt8}) where {N, T}
    return Category{T}(StaticString(vec))
end

function name(category::Category)
    category.str
end

function ==(lhs::Category, rhs::Category)
    name(lhs) == name(rhs)
end

function Base.length(category::Category)
    return length(name(category))
end

struct SerializedCategory
    chars::Vector{UInt8}
    bounds::Vector{Int64}
end

# serialization for data storage

function serialize(cats::Vector{Category{T}}) where {T}
    len = sum(length(e) for e in cats)
    chars = Vector{UInt8}(undef, len)
    bounds = Vector{Int64}(undef, length(cats))
    nchar = 1
    for (i, e) in enumerate(cats)
        bounds[i] = nchar
        for c in name(e)
            chars[nchar] = c
            nchar += 1
        end
    end

    return SerializedCategory(chars, bounds)
end

function deserialize(scats::SerializedCategory)
    chars = scats.chars
    bounds = scats.bounds
    cats = [Category{Element}(StaticString("")) for i in 1:length(bounds)]
    for i in 1:length(bounds)-1
        s, f = bounds[i], bounds[i+1]-1
        cname = Tuple(c for c in view(chars, s:f))
        cats[i] = Category{Element}(cname)
    end
    cname = Tuple(c for c in view(chars, bounds[end-1]:length(chars)))
    cats[end] = Category{Element}(cname)

    return cats
end
