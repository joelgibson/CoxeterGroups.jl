# CoxeterGroups.jl Documentation

Coxeter groups and their elements are subtypes of the abstract types [`CoxGrp`](@ref) and [`CoxElt`](@ref).
Currently there are two (incomplete) implementations, which can be constructed using either [`CoxeterGroup`](@ref) or [`coxeter_group_min`](@ref).
The constructor takes a matrix, which can be either a generalised Cartan matrix, or a Coxeter matrix (currently `coxeter_group_min` cannot accept a Coxeter matrix), and returns a group object and list of generators.
For example:

```julia
W, (s, t) = coxeter_group_min([2 -1; -1 2])
s*t*s == t*s*t
```

## Creating a group

There are various implementations of Coxeter groups, each represented by a concrete type which is a subtype of the abstract [`CoxGrp`](@ref).
Use one of the constructors below to create a group, depending on which implementation is desired.
For example:

```julia
# A Coxeter group of type A2, implemented using minimal roots.
W, (s, t) = coxter_group_min([2, -1; -1, 2])

# The symmetric group S3, implemented via permutations.
S, (u, v) = symmetric_group(3)
```

Here is a summary of the creation functions - there is more information on specific implementations in the [Implementations](@ref) section.
```@docs
coxeter_group_min
symmetric_group
```


## Operations on the group

A Coxeter group is a subtype of the abstract [`CoxGrp`](@ref) type.
Coxeter groups compare equal by reference.
The following operations are supported on all group implementations:

```@docs
CoxeterGroups.rank(::CoxGrpMin)
CoxeterGroups.generators(::CoxGrpMin)
Base.one(::CoxGrpMin)
CoxeterGroups.is_finite(::CoxGrpMin)
CoxeterGroups.longest_element(::CoxGrp)
```

## Operations on group elements

A Coxeter group element is a subtype of the abstract [`CoxElt`](@ref) type.
Coxeter group elements are immutable value types: they compare equal if they represent the same group element, they can be used as keys in hash tables, and so on.
The following operations are supported on all Coxeter group element implementations:

```@docs
Base.parent(::CoxEltMin)
Base.isone(::CoxEltMin)
Base.:(*)(::CoxEltMin, ::CoxEltMin)
Base.inv(::CoxEltMin)
Base.:(^)(::CoxElt, ::Integer)
CoxeterGroups.length(::CoxEltMin)
Base.sign(::CoxElt)
CoxeterGroups.is_right_descent(::CoxEltMin, ::Integer)
CoxeterGroups.is_left_descent(::Integer, ::CoxEltMin)
CoxeterGroups.right_multiply(::CoxEltMin, ::Integer)
CoxeterGroups.left_multiply(::Integer, ::CoxEltMin)
CoxeterGroups.short_lex(::CoxEltMin)
CoxeterGroups.inverse_short_lex(::CoxEltMin)
```

## Internals

A Coxeter group implementation consists of a new pair of concrete types subtyping [`CoxGrp`](@ref) and [`CoxElt`](@ref) respectively.

```@docs
CoxGrp
CoxElt
```

When adding a new concrete implementation, a small set of functions must be implemented explicitly, then generic implementations of other functions will begin working.
The concrete implementation should specialise these generic implementations where it is possible to take advantage of specific internal structure.

- The group type must implement a constructor such as those in [Creating a group](@ref), and the functions [`coxeter_matrix`](@ref) [`rank`](@ref), [`one`](@ref), [`generators`](@ref), and [`is_finite`](@ref).
- The element type must implement [`parent`](@ref), `==`, `hash`, [`is_left_descent`](@ref), [`is_right_descent`](@ref), [`left_multiply`](@ref), and [`right_multiply`](@ref).

The following functions will then be defined:

- For the group type: [`longest_element`](@ref).
- For the element type: [`isone`](@ref), [`length`](@ref), [`sign`](@ref), [`short_lex`](@ref), [`inverse_short_lex`](@ref), [`*`](@ref), [`inv`](@ref), [`^`](@ref).



## Implementations

### Minimal roots implementation

The _minimal roots_ implementation of a Coxeter group represents group elements by words in the Coxeter generators in ShortLex normal form.
Operations on these words are done with the assistance of a _minimal root reflection table_, which needs to be constructed when the group is created.

Currently this implementation can handle any Coxeter group defined by a generalised Cartan matrix, but this is only because generating the minimal root reflection table is more difficult for general Coxeter groups.
The same underlying algorithms will work once given a new table.

Performance characteristics for a Coxeter group of rank ``r``:

- The group itself takes ``O(r M)`` space, where ``M`` is the number of minimal roots.
  For finite types, the minimal roots are exactly the positive roots, and for affine types, there are twice as many minimal roots as there are positive roots for the corresponding finite system.
  The number of roots in a finite system tends to be ``O(r^2)`` in the rank ``r``.
- Length is an ``O(1)`` operation.
- ShortLex normal form should be ``O(1)``, but is ``O(l(w))`` since it returns a defensive copy of the internal word.
- Testing a left or right descent is ``O(l(w))``.
- Multiplication on the right ``xs`` by a simple generator ``s`` is ``O(l(x))``.
- Multiplication on the left ``sx`` by a simple generator ``s`` is ``O(l(x)^2)``.
- Multiplication ``xy`` is ``O(l(y)(l(x) + l(y)))``.
- Inversion ``x^{-1}`` is ``O(l(x)^2)``.

```@docs
CoxGrpMin
CoxeterGroups.create_reflection_table_gcm
```

### Symmetric groups

The symmetric group ``S_n`` on ``n`` letters is the group of permutations of the set ``[n] = \{1, \ldots, n\}``, whose standard Coxeter generators are ``s_1, \ldots, s_{n-1}`` where ``s_i = (i, i+1)`` is an adjacent transposition.
The symmetric group may be constructed using the [`symmetric_group`](@ref) function.


A permutation ``\pi \colon [n] \to [n]`` is represented as a pair of arrays recording both ``\pi`` and its inverse ``\pi^{-1}``:
```math
    \pi \mapsto (\mathtt{fwd} = [\pi(1), \ldots, \pi(n)], \mathtt{inv} = [\pi^{-1}(1), \ldots, \pi^{-1}(n)]).
```

Complexity of operations in terms of ``n``:

- Creation of the group is ``O(1)``.
- Memory usage for each element is ``O(n)``.
- Testing a left or right descent is ``O(1)``.
- Multiplication is ``O(n)``.
- Calculating the sign of a permutation is ``O(n)``
- Taking Coxeter length takes ``O(n^2)`` time (this could be improved to ``O(n \log n)`` by a standard divide-and-conquer trick to count inversions).
