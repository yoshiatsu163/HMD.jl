#####
##### StaticString definition
#####

struct StaticString
    #str::NTuple{N, UInt8}
    str::Vector{UInt8}
end

function StaticString(str::AbstractString)
    #sstr = NTuple{N, UInt8}(UInt8(c) for c in str)

    return StaticString([UInt8(c) for c in str])
end

#function StaticString(vec::AbstractVector{Char})
#    return StaticString{length(vec)}(Tuple(vec))
#end

function string(sstr::StaticString)
    return UInt8.(sstr.str) |> join
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

#function Category{T}(vec::Vector{Char}) where {T}
#    return Category{T}(vec)
#end

function Category{T}(vec::AbstractVector{UInt8}) where {T}
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

# serialization for data storage

"""
    serialize(strings)

convert Vector{String} to chars::Vector{UInt8} and bounds::Vector{Int64}.
The first vector contains all the characters in the strings,
and the second vector contains the starting index of each string.

"""
function serialize(strings::Vector{String})
    chars = Vector{UInt8}(undef, sum(length(str) for str in strings))
    bounds = Vector{Int64}(undef, length(strings))
    nchar = 1
    for i in 1:length(strings)
        bounds[i] = nchar
        for c in codeunits(strings[i])
            chars[nchar] = c
            nchar += 1
        end
    end
    push!(bounds, nchar)

    return chars, bounds
end

"""
    deserialize(chars, bounds)

convert chars::Vector{UInt8} and bounds::Vector{Int64} to Vector{string}.
`char` contains all the characters in the strings,
and `bounds` contains the starting index of each string in `chars`.
"""
function deserialize(chars::Vector{UInt8}, bounds::Vector{Int64})
    return map(1:length(bounds)-1) do i
        String(chars[bounds[i]:bounds[i+1]-1])
    end
end
