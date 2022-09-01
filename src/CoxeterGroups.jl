# By T. Schmit (2021)

module CoxeterGroups

export CoxeterGroup, CoxeterElement

# Abstract supertypes
export CoxGrp, CoxElt

# Group functions
export coxeter_matrix, rank, generators, longest_element, is_finite

# Element functions
export length, sign, is_right_descent, is_left_descent, right_multiply, left_multiply, short_lex, inverse_short_lex


"""The abstract supertype for Coxeter groups."""
abstract type CoxGrp end

"""The abstract supertype for Coxeter group elements."""
abstract type CoxElt end

include("CoxeterGroupData.jl")
include("CoxMin.jl")
include("SymmetricGroup.jl")

# Default implementations of some functions.

"""
    isone(w)

Return true if the group element is the identity."""
Base.isone(w::CoxElt) = w == one(parent(w))

"Return the numerically least left descent of w, or 0 if there is no left descent (iff w = id)."
function first_left_descent(w::CoxElt)
    for i in 1:rank(parent(w))
        if is_left_descent(i, w)
            return i
        end
    end
    return 0
end

"Return the numerically least right descent of w, or 0 if there is no right descent (iff w = id)."
function first_right_descent(w::CoxElt)
    for i in 1:rank(parent(w))
        if is_right_descent(w, i)
            return i
        end
    end
    return 0
end

"""
    length(w::CoxElt)

The Coxeter length of ``w``: the length of a shortest expression in the Coxeter generators for ``w``."""
function Base.length(w::CoxElt)
    length = 0
    while true
        s = first_left_descent(w)
        s == 0 && return length
        length += 1
        w = left_multiply(s, w)
    end
end

"""
    sign(w::CoxElt) -> Int

The sign homomorphism ``(-1)^{l(w)}``.
"""
Base.sign(w::CoxElt) = length(w) % 2 == 0 ? 1 : -1

"""
    short_lex(w)

Return the ShortLex normal form of w: the lexicographically least reduced expression for w.
"""
function short_lex(w::CoxElt)
    word = []
    while true
        s = first_left_descent(w)
        s == 0 && return word
        push!(word, s)
        w = left_multiply(s, w)
    end
end

"""
    inverse_short_lex(w)

Return the InverseShortLex normal form of `w`: the lexicographically least reduced expression for ``w``
when read backwards, or equivalently the reverse of the ShortLex normal form of ``w^{-1}``.
"""
function inverse_short_lex(w::CoxElt)
    word = []
    while true
        s = first_right_descent(w)
        s == 0 && return reverse!(word)
        push!(word, s)
        w = right_multiply(w, s)
    end
end

"""
    *(x, y)

Take the product of two Coxeter group elements."""
function Base.:(*)(x::CoxElt, y::CoxElt)
    parent(x) == parent(y) || error("Group elements $x and $y do not come from the same parent.")
    while true
        s = first_left_descent(y)
        s == 0 && return x
        x = right_multiply(x, s)
        y = left_multiply(s, y)
    end
end

# The product over one element sometimes comes up, eg when doing *(gens...) to get a Coxeter element.
Base.:(*)(x::CoxElt) = x

"""
    inv(x)

The inverse of a Coxeter group element."""
function Base.inv(x::CoxElt)
    inv = one(parent(x))
    while true
        s = first_left_descent(x)
        s == 0 && return inv
        x = left_multiply(s, x)
        inv = left_multiply(s, inv)
    end
end

"""
    ^(x::CoxElt, e::T) where T <: Integer

Return the power ``x^e``."""
function Base.:(^)(x::CoxElt, e::T) where T <: Integer
    e == 0 && return one(parent(x))
    if e < 0
        x = inv(x)
        e = -e
    end

    e == 1 && return x

    # Exponentiate by squaring
    acc = one(parent(x))
    pow2 = x
    while e != 0
        if e % 2 == 1
            acc *= pow2
        end
        e = div(e, 2)
        pow2 *= pow2
    end
    return acc
end

"""
    longest_element(W::CoxGrp)

Return the longest element if ``W`` is finite, otherwise error."""
function longest_element(grp::CoxGrp)
    is_finite(grp) || error("longest_element only works for finite Coxeter groups")
    w = one(grp)
    while true
        finished = true
        for i in 1:rank(grp)
            if !is_right_descent(w, s)
                w = right_multiply(w, s)
                finished = false
            end
        end

        finished && return w
    end
end


"""
    CoxeterGroup

A **Coxeter group** is a group **G** together with a generating subset **S**, whose elements are subject to relations:
```math
\\begin{aligned}
s² &= 1 \\quad \\text{ for } s∈S \\\\

(st)^{m_{s,t}} = 1 \\quad \\text{ for } s≠t∈S,\\; m_{s,t}∈\\mathbb{N}_{≥2} \\\\
\\end{aligned}
```
A CoxeterGroup can be defined by its Coxeter matrix M=(m_{s,t})_{s,t∈S} and the generating subset **S**.
Elements in CoxeterGroup are represented by their expression in InverseShortLex form.
Once the CoxeterGroup has been established, you can calculate on it using "*".

# Example
```julia-repl
julia> G,(a,b,c) = CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"])
(CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"]), CoxeterElement[[1], [2], [3]])
julia> w = a*c
c*a
julia> w*c
a
```
"""
struct CoxeterGroup <: CoxGrp
    M::Matrix{Int64}
    S::Vector{String}

    function CoxeterGroup(M::Matrix{Int64})
        if is_gcm(M)
            M = gcm_to_coxeter_matrix(M)
        end
        is_coxeter_matrix(M) || throw(ArgumentError("M has to be a Coxeter matrix"))

        S = ["<$i>" for i=1:size(M,1)]
        CG = new(M,S)
        return CG, [CoxeterElement(CG, [i], true) for i=1:length(S)]
    end

    function CoxeterGroup(M::Matrix{Int64}, S::Array{String, 1})
        if is_gcm(M)
            M = gcm_to_coxeter_matrix(M)
        end
        is_coxeter_matrix(M) || throw(ArgumentError("M has to be a Coxeter matrix"))
        length(S) == size(M,1) || throw(ArgumentError("For M of size $(size(M)), S should have length $(size(M,1))"))

        CG = new(M,S)
        return CG, [CoxeterElement(CG, [i], true) for i=1:length(S)]
    end
end


"""
    CoxeterElement

A **Coxeter element** is an element of a Coxeter group. We use CoxeterElement in order to facilitate the calculation in a given Coxeter group.
A CoxeterElement is always represented by its InverseShortLex form.
"""
struct CoxeterElement <: CoxElt
    G::CoxeterGroup
    w::Array{Int,1}
    """
        CoxeterElement(G::CoxeterGroup, w::Array{Int,1}, isNormal=false::Bool)

    Returns the CoxeterElement corresponding to S[w] for S the generating set of G.
    If w is already in InverseShortLex form, you can set isNormal=true in order to save on computation.

    # Example
    ```julia-repl
    julia> G,S = CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"])
    (CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"]), CoxeterElement[[1], [2], [3]])
    julia> w = CoxeterElement(G, [1,2,3])
    abc
    julia> v = CoxeterElement(G, [3,2,3])
    bcb
    ```
    """
    function CoxeterElement(G::CoxeterGroup, w::Array{Int,1}, isNormal=false::Bool)
        if isNormal
            new(G, w)
        else
            new(G, NF(w,G))
        end
    end
end

#This function is only for convenience. It can only be called to produce the identity element
function CoxeterElement(G::CoxeterGroup, w::Array{Any,1}, isNormal=false::Bool)
    isempty(w) || throw(ArgumentError("w has to be of type Array{Int,1}; I received an element of type $(typeof(w))"))
    return CoxeterElement(G, Int[], true)
end


function Base.show(io::IO, mime::MIME"text/plain", w::CoxeterElement)
    output = ""
    for s in w
        output *= w.G.S[s]
    end
    if length(w) > 0
        print(io, output)
    else
        print(io, "<>")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", C::CoxeterGroup)
    println(io, "Coxeter Group with Coxeter Matrix:")
    show(io, mime, C.M)
    println(io, "\nand generating Set:")
    println(io, C.S)
end

"""The identity element of the Coxeter group."""
Base.one(C::CoxeterGroup) = CoxeterElement(C, Int[], true)

"""The Coxeter generators."""
generators(C::CoxeterGroup) = [CoxeterElement(C, [i], true) for i=1:rank(C)]


Base.length(w::CoxeterElement) = length(w.w)

Base.getindex(w::CoxeterElement, i) = getindex(w.w, i)


function Base.insert!(w::CoxeterElement, i::Int, s::Int)
    insert!(w.w, i, s)
end

function Base.deleteat!(w::CoxeterElement, i::Int)
    deleteat!(w.w, i)
end

Base.copy(w::CoxeterElement) = CoxeterElement(w.G, copy(w.w), true)

Base.:(==)(w::CoxeterElement, x::CoxeterElement) = w.G == x.G && w.w == x.w

Base.hash(w::CoxeterElement, h::UInt) = hash(w.w, h)

"""
    <(x::CoxeterElement, y::CoxeterElement)

Returns true if x<y according to the InverseShortLex ordering.
"""
function Base.:(<)(x::CoxeterElement, y::CoxeterElement)
    x.G.M==y.G.M || throw(ArgumentError("x and y do not belong to the same Coxeter group"))
    if length(x) < length(y)
        return true
    elseif length(x) > length(y)
        return false
    else
        i = length(x)
        while i>0
            if x[i] == y[i]
                i-=1
            else
                return (x[i]<y[i])
            end
        end
        return false
    end
end

"""
    >(x::CoxeterElement, y::CoxeterElement)

Returns true if x>y according to the InverseShortLex ordering.
"""
function Base.:(>)(x::CoxeterElement, y::CoxeterElement)
    return <(y,x)
end

"""
    *(x::CoxeterElement, y::CoxeterElement)

The operation of the Coxeter group. Returns the CoxeterElement in InverseShortLex form of x⋅y.
"""
function Base.:(*)(x::CoxeterElement, y::CoxeterElement)
    x.G.M==y.G.M || throw(ArgumentError("x and y do not belong to the same Coxeter group"))
    return NF(x,y)
end

"""
    coxeter_matrix(CG::CoxeterGroup)

Returns the associated Coxeter matrix for CG
"""
coxeter_matrix(CG::CoxeterGroup) = deepcopy(CG.M)

"""
    parent(w::CoxeterElement)

Returns the Coxeter group in which w currently resides.
"""
Base.parent(w::CoxeterElement) = w.G

"""
    rank(CG::CoxeterGroup)

Returns the size of the generating set for CG.
"""
rank(CG::CoxeterGroup) = length(CG.S)




#returns the normal form of S[x] in CG
function NF(x::Array{Int,1}, CG::CoxeterGroup)
    if length(x)>1
        return NF(x[1:end-1], CoxeterElement(CG, [x[end]], true))
    else
        return CoxeterElement(CG, x, true)
    end
end

#returns the normal form of x*y
function NF(x::CoxeterElement, y::CoxeterElement)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end

#returns the normal form of x*y
function NF(x::Array{Int,1}, y::CoxeterElement)
    result = copy(y)
    for i = length(x):-1:1
        result = NF!(x[i], result)
    end
    return result
end


function NF(s::Int, w::CoxeterElement)
    return NF!(s, copy(w))
end


function NF!(s::Int, w::CoxeterElement)
    n = length(w)
    l = 1
    r = 1
    σ = s
    while r <= n
        c = w[r]
        r = r + 1
        if c >= σ
            σ = c
            t = Exchange(s, CoxeterElement(w.G, w.w[l:r-1], true))
            if t != nothing
                if s == w[l]
                    deleteat!(w,l)
                    return(w)
                else
                    s = t
                    l = r
                    σ = t
                end
            end
        end
    end
    insert!(w, l, s)
    return(w)
end


function Exchange(s::Int, w::CoxeterElement)
    CG = w.G
    M = CG.M
    n = length(w)
    if n == 1
        if w[1] == s
            return s
        elseif w[1]<s || M[w[1],s]!=2
            return nothing
        else
            return s
        end
    end
    if n == M[w[1], s]-1
        i = 2
        while i<=n && (w[i]==s || w[i]==w[1])
            i += 1
        end
        if i == n + 1
            if w[n] > w[n-1]
                return w[n-1]
            else
                return nothing
            end
        end
    else
        sy = CoxeterElement(CG, pushfirst!(w[1:n-1], s), true)
        z = NF(w[1], sy) #s_1sy
        if length(z) < n
            m = 1
            while m<length(z) && z[m+1] == sy[m+1]
                m = m + 1
            end
            ω = NF(vcat(w[1], s, w[1:m-1]), CoxeterElement(CG, w[m+1:n], true))

            for t in ω.w
                if t > ω.w[end]
                    return ω.w[end]
                end
            end
        end
        return nothing
    end
end

end # module
