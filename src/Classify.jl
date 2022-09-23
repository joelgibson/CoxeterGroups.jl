using Match

export classify_coxeter_matrix, classify_gcm, is_finite_type, cartan_matrix, root_system, is_finite_type, positive_roots, positive_coroots, coxeter_system, degrees, order, coxeter_number, number_of_reflections, exponents

#=
    Classification of finite and affine-type Dynkin diagrams, and finite-type Coxeter systems.

The following table is a list of all finite and affine type Dynkin diagrams and finite type Coxeter systems. The I(GCM)
and I(Cox) columns specify which types should be taken to make an irredundant list of finite and affine type GCMs or
Coxeter systems - when a Coxeter type should be omitted, an isomorphic type is placed in there instead. The Shape column
classifies the underlying undirected graph, which is either a path, cycle, or tree. When the underlying graph is a tree,
the vertex degrees which are ≥ 3 are listed. The Bond column is the multiset of Coxeter bonds, excluding m_st = 2 or 3.
The ~ symbol means affinisation, and the @ symbol means dual affinisation.

Type         I(GCM)    I(Cox)   Shape       Bond
----         ------    ------   -----       ----
An           n ≥ 1     n ≥ 1    Path        -
Bn           n ≥ 3     Cn       Path        4
Cn           n ≥ 2     n ≥ 2    Path        4
Dn           n ≥ 4     n ≥ 4    Tree(3)     -
E6           Yes       Yes      Tree(3)     -
E7           Yes       Yes      Tree(3)     -
E8           Yes       Yes      Tree(3)     -
F4           Yes       Yes      Path        4
G2           Yes       Yes      Path        6

Hn           n=2,3,4   n=2,3,4  Path        5
I2(m)        m ≥ 7     m ≥ 7    Path        m

A~1          Yes       Yes      Path        ∞
A~n          n ≥ 2     n ≥ 2    Cycle       -
B~n          n ≥ 3     n ≥ 3    Tree(3)     4
C~n          n ≥ 2     n ≥ 2    Path        4,4                             Bonds point in.
D~4          Yes       Yes      Tree(4)     -
D~n          n ≥ 5     n ≥ 5    Tree(3,3)   -
E~6          Yes       Yes      Tree(3)     -
E~7          Yes       Yes      Tree(3)     -
E~8          Yes       Yes      Tree(3)     -
F~4          Yes       Yes      Path        4
G~2          Yes       Yes      Path        6

BC@1         Yes       A~1      Path        ∞       Kac: A2(2)
BC@n         n ≥ 2     C~n      Path        4,4     Kac: A_{2(n - 1)}(2)    Bonds point same direction.
B@n          n ≥ 3     B~n      Tree(3)     4       Kac: A_{2n - 1}(2)
C@n          n ≥ 2     C~n      Path        4,4     Kac: D_{n-1}(2)         Bonds point out.
F@4          Yes       F~4      Path        4,4     Kac: E6(2)
G@2          Yes       G~2      Path        6       Kac: D4(3)

Types are named as follows:
- Finite type GCMs are a pair (type, rank), for example (:C, 4)
- Affine type GCMs are a triple (type, rank, kind), where (type, rank) specifies a finite type irreducible GCM, and kind
  is either :aff (the usual affinisation inside Aff(coroot space)) or :dualaff (affinisation inside Aff(root space)).
- The other finite Coxeter groups are named (:H, n) for n in 2, 3, 4, and (:I, 2, m).
- Something not on the classification list is classified as (:Unknown, rank).

The differences between types Bn and Cn, and between :aff and :dualaff, are only visible on the Cartan matrix rather
than the Coxeter matrix. When classifying Coxeter matrices, type :C will be preferred to :B, and :aff to :dualaff.
=#

"""
    classify_coxeter_matrix(mat)

Classify a Coxeter matrix as one of the finite or affine type Coxeter systems. A list of pairs is returned, with each
pair giving a classification such as (:A, 4) or (:D, 6, :aff), together with an ordered list of vertices making up that
component in the Kac convention.
"""
function classify_coxeter_matrix(coxmat::AbstractMatrix{T}) where {T <: Integer}
    is_coxeter_matrix(coxmat) || error("Argument was not a Coxeter matrix")

    # Rank of the large Coxeter matrix.
    coxrank = size(coxmat)[1]

    # Adjacency list and degrees for the underlying undirected graph.
    adj = Dict(s => Int[t for t in 1:coxrank if s != t && coxmat[s, t] != 2] for s in 1:coxrank)
    degrees = [length(adj[s]) for s in 1:coxrank]


    # Create a list of the connected components of the underlying undirected graph.
    components = Vector{Vector{Int64}}()
    seen = zeros(Bool, coxrank)
    for start in 1:coxrank
        seen[start] && continue

        component = [start]
        seen[start] = true
        pos = 1
        while pos <= length(component)
            node = component[pos]
            pos += 1
            for t in 1:coxrank
                if !seen[t] && coxmat[node, t] != 2 && coxmat[node, t] != 1
                    push!(component, t)
                    seen[t] = true
                end
            end
        end
        push!(components, component)
    end

    # Perform a breadth-first search from a point, returning a triple (vertices, dist, pred) of the vertices visited
    # in BFS order, their distances from the source, and a dictionary mapping each vertex to its predecessor on a
    # shortest path to the source.
    function bfs(start)
        dist = Dict(start => 0)
        pred = Dict(start => start)
        order = [start]
        pos = 1
        while pos <= length(order)
            s = order[pos]
            pos += 1
            for t in 1:coxrank
                if s != t && coxmat[s, t] != 2 && !haskey(dist, t)
                    dist[t] = dist[s] + 1
                    pred[t] = s
                    push!(order, t)
                end
            end
        end
        return (order, dist, pred)
    end

    # Given a pred dictionary output by BFS, and a starting vertex v, return the shortest path [v, ..., BFS source.]
    function getpath(pred, v)
        path = [v]
        while true
            v = path[end]
            pred[v] == v && return path
            push!(path, pred[v])
        end
    end

    # Given a vertex in a cycle, return a list of vertices in that cycle, starting from start, in cycle order.
    function getcycle(start)
        degrees[start] == 2 || error("Degree of start should be 2 for getcycle")
        order = [start]

        prev = start
        next = minimum(adj[start])
        while next != start
            push!(order, next)
            degrees[next] == 2 || error("Not a cycle")
            (prev, next) = (next, [s for s in adj[next] if s != prev][1])
        end

        return order
    end


    # Classify a component, returning something like (:A, [4, 1]) or (:Unknown, []).
    function classify_component(comp)
        # The rank of this irreducible component.
        rank = length(comp)

        # List, with multiplicity, the bonds m ∈ {0, 4, 5, ...}
        multbonds = sort(Int[m for m in [coxmat[comp[s], comp[t]] for s in 1:rank for t in s+1:rank] if m != 2 && m != 3])

        # A connected graph is a tree iff |E| = |V| - 1, and a cycle iff deg(v) = 2 for all v ∈ V.
        nedges = sum(degrees[s] for s in comp) ÷ 2
        is_tree = rank - 1 == nedges
        is_cycle = all(degrees[s] == 2 for s in comp)

        # The only non-tree graph we accept is A~n for n ≥ 2, which is a cycle.
        if !is_tree
            # If we have a cycle where every edge has bond multiplicity m = 3, then we have affine type A.
            is_cycle && multbonds == [] && return ((:A, rank-1, :aff), getcycle(minimum(comp)))

            # Otherwise, unknown type.
            return ((:Unknown, rank), comp)
        end

        # If the tree is made up of one vertex, then we have type A1.
        nedges == 0 && return ((:A, 1), comp)


        # List the leaves of the tree, and for each leaf, the bond multiplicity incident on it.
        leaves = sort([s for s in comp if degrees[s] == 1])
        leaf_bond = Dict(s => [coxmat[s, t] for t in comp if s != t && coxmat[s, t] != 2][1] for s in leaves)

        # The degrees of the component which are ≥ 3.
        bigdegrees = sort(Int[degrees[s] for s in comp if degrees[s] >= 3])

        # A path is a tree where all degrees are ≤ 2, so bigdegrees will be empty.
        if bigdegrees == []
            # Reorder the leaves so that the one with the lower bond multiplicity comes first.
            sort!(leaves; by=leaf -> leaf_bond[leaf])

            #                       mleft               mright
            # Let the path be left ------- (some path) -------- right, where mleft ≤ mright. Reorder comp so that
            # vertices are in path order [left, ..., right].
            (left, right) = leaves
            (comp, _, _) = bfs(left)
            (comp[1] == left && comp[rank] == right) || error("Comp should go from left to right now.")

            # The bond multiplicity incident on the left and right leaves.
            mleft = leaf_bond[left]
            mright = leaf_bond[right]


            ## Finite-type Coxeter systems which are paths.

            # An: No bonds of multiplicity ≂̸ 3.
            multbonds == [] && return ((:A, rank), comp)

            # Cn: Extra bond multiplicities {4}, and that bond is incident on the right vertex.
            multbonds == [4] && mright == 4 && return ((:C, rank), comp)

            # F4: Extra bond multiplicities {4}, rank 4, and both leaves have simple bonds.
            multbonds == [4] && rank == 4 && mleft == mright == 3 && return ((:F, 4), comp)

            # G2: Extra bond multiplicities {6}, rank 2.
            multbonds == [6] && rank == 2 && return ((:G, 2), comp)

            # Hn: Extra bond multiplicities {5}, that 5 is incident on the right vertex, and rank is 2, 3, or 4.
            multbonds == [5] && mright == 5 && 2 <= rank <= 4 && return ((:H, rank), comp)

            # I2(m): Extra bond multiplicities {m} for some m ≥ 7, rank 2.
            m = (length(multbonds) == 1) ? multbonds[1] : -1
            rank == 2 && m >= 7 && return ((:I, 2, m), comp)


            ## Affine-type Coxeter systems which are paths.

            # A~1: Extra bond multiplicities {∞}, rank 2.
            rank == 2 && multbonds == [0] && return ((:A, 1, :aff), comp)

            # C~n: Extra bond multiplicities {4, 4}, each incident on a leaf.
            # Coxeter graph: (n+1) == 1 -- 2 ... -- (n-2) -- (n-1) == n
            multbonds == [4, 4] && mright == mleft == 4 && return ((:C, rank-1, :aff), [comp[2:n+1]; comp[1]])

            # F~4: Extra bond multiplicities {4}, not incident on a leaf, and rank 5. In Kac, the bond in F4 points to
            # the right, and so in the affinisation we should have [3, 3, 4, 3] as the path bonds.
            if multbonds == [4] && mleft == 3 && mright == 3
                # Kac labelling should have 5 -- 1 -- 2 == 3 -- 4
                comp = coxmat[comp[3], comp[4]] == 4 ? comp : reverse(comp)
                return ((:F, 4, :aff), [comp[2:5]; comp[1]])
            end

            # G~2: Extra bond multiplicities {6}, incident on the right leaf, with rank 3.
            # Currently our path looks like 1 -- 2 ≡≡ 3, we need to return 3 -- 1 ≡≡ 2
            multbonds == [6] && rank == 3 && mright == 6 && return ((:G, 2, :aff), [comp[2:3]; comp[1]])

            ## No more finite or affine-type Coxeter systems which are paths.
            return ((:Unknown, rank), comp)
        end

        # Now for those trees with a unique vertex of degree 3, and no other vertices of higher degree.
        if bigdegrees == [3]
            # A tree with a unique vertex of degree 3 is essentially three paths joined to a special "star" point.
            # It has three leaves: let these leaves be p, q, r ordered by distance from the star point, with ties
            # broken according to their incident bond.

            # Perform a BFS to get distances from the star vertex.
            star = [s for s in comp if degrees[s] == 3][1]
            (_, dist, pred) = bfs(star)

            # Reorder the leaves by the criterion above, and take vectors of the resulting distances and bond mults.
            sort!(leaves; by=leaf -> (dist[leaf], leaf_bond[leaf]))
            leaf_dists = [dist[leaf] for leaf in leaves]
            leaf_bonds = [leaf_bond[leaf] for leaf in leaves]

            ## Finite-type Coxeter systems which are trees.

            # Dn: All simple edges, with at least two leaves having distance 1 to the star vertex.
            if multbonds == [] && leaf_dists[1:2] == [1, 1]
                # In type D4 in Kac, the star vertex is the second vertex.
                rank == 4 && return ((:D, 4), [star ; leaves])

                # Otherwise, BFS from that unique leaf at distance ≥ 2 from the star vertex to get the ordering.
                (comp, _, _) = bfs(leaves[3])
                return ((:D, rank), comp)
            end

            # En: All simple edges, ranks 6, 7, 8, with leaves at distances [1, 2, ?] from the star vertex.
            if multbonds == [] && 6 <= rank <= 8 && leaf_dists[1:2] == [1, 2]
                # For E6 and E7, we want a path 1 -- 2 -- 3star -- 4 -- 5 ( -- 6 ), with the star leaf attached to 3.
                # For E8 however, we want 1 -- 2 -- 3 -- 4 -- 5star -- 6 -- 7 to fit with Kac' labelling.
                order = (
                    (rank == 6 || rank == 7)
                    ? [getpath(pred, leaves[2]) ; getpath(pred, leaves[3])[end-1:-1:1] ; leaves[1]]
                    : [getpath(pred, leaves[3]) ; getpath(pred, leaves[2])[end-1:-1:1] ; leaves[1]]
                )
                return ((:E, rank), order)
            end


            ## Affine-type Coxeter systems which are trees.

            # B~n: Extra bond multiplicities {4}, two leaves incident on simple bonds at distance 1, and the last
            # leaf incident on the bond with multiplicity 4.
            multbonds == [4] && leaf_dists[1:2] == [1, 1] && leaf_bonds[3] == 4 && (
                return ((:B, rank - 1, :aff), [leaves[1] ; getpath(pred, leaves[3])[end:-1:1] ; leaves[2]])
            )

            # E~6: All bonds m = 3, and each leaf at distance 2.
            if multbonds == [] && leaf_dists == [2, 2, 2]
                # Extended diagram is
                #             aff
                #              |
                #              6
                #              |
                #  1 --- 2 --- 3 --- 4 --- 5
                nodes = [
                    getpath(pred, leaves[1]);
                    getpath(pred, leaves[2])[end:-1:2];
                    getpath(pred, leaves[3])[end:-1:2]
                ]
                return ((:E, 6, :aff), nodes)
            end

            # E~7: All bonds m=3, leaf distances 1, 3, and 3.
            if multbonds == [] && leaf_dists == [1, 3, 3]
                # Extended diagram is
                #                     7
                #                     |
                # aff --- 1 --- 2 --- 3 --- 4 --- 5 --- 6
                (aff, left...) = getpath(pred, leaves[2])
                right = getpath(pred, leaves[3])[end:-1:2]
                return ((:E, 7, :aff), [left ; right ; leaves[1] ; aff])
            end

            # E~8: All bonds m=3, leaf distances 1, 2, 5.
            if multbonds == [] && leaf_dists == [1, 2, 5]
                # Extended diagram is
                #                                 8
                #                                 |
                # aff --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 --- 7
                (aff, left...) = getpath(pred, leaves[3])
                right = getpath(pred, leaves[2])[end:-1:2]
                return ((:E, 8, :aff), [left ; right ; leaves[1] ; aff])
            end


            ## No more finite or affine type Dynkin diagrams of this kind.
            return ((:Unknown, rank), comp)
        end

        # D~4: Branching degrees {4}, simple bonds throughout, rank 5.
        if bigdegrees == [4] && multbonds == [] && rank == 5
            # The 4-valent vertex comes second.
            star = [s for s in comp if degrees[s] == 4][1]
            return ((:D, 4, :aff), [leaves[1] ; star ; leaves[2:end]])
        end

        # D~n for n ≥ 5: Branching degrees {3, 3}, simple bonds, each leaf incident on a branching vertex.
        leaf_neighbour_degrees = [degrees[adj[leaf][1]] for leaf in leaves]
        if bigdegrees == [3, 3] && multbonds == [] && leaf_neighbour_degrees == [3, 3, 3, 3]
            # Extended diagram is
            #       aff             n
            #        |              |
            #  1 --- 2 --- ... --- n-2 --- n-1
            (order, _, _) = bfs(leaves[1])
            affleaf = [leaf for leaf in leaves if leaf in adj[leaves[1]]]
            return ((:D, rank - 1, :aff), [[s for s in order if s != affleaf] ; affleaf])
        end

        # This exhausts all finite and affine type Coxeter systems.
        return ((:Unknown, rank), comp)
    end

    return [classify_component(component) for component in components]
end


function classify_gcm(gcm::AbstractMatrix{T}) where {T <: Integer}
    is_gcm(gcm) || error("Was not given a GCM")

    # First run the classification of the underlying Coxeter graph, then check bond directions to distinguish
    # between dual cases. This relies on the convention we know which is being returned in the classification: eg if
    # a Coxeter graph is classified as type C, then the double bond is at the end.
    cox_components = classify_coxeter_matrix(gcm_to_coxeter_matrix(gcm))

    convert_component(type, comp) = @match type begin
        # Cases where we can just forward the classification of the Coxeter matrix.
        (:Unknown, n)                   => (type, comp)
        (:A, n)                         => (type, comp)
        (:A, n, :aff), if n >= 2 end    => (type, comp)
        (:D, n)                         => (type, comp)
        (:D, n, :aff)                   => (type, comp)
        (:E, n)                         => (type, comp)
        (:E, n, :aff)                   => (type, comp)

        # Non-simply laced finite types. Note the Coxeter classifier will never return type B, only type C.
        (:C, n) => gcm[comp[n], comp[n-1]] == -2 ? ((:B, n), comp) : ((:C, n), comp)
        (:F, 4) => gcm[comp[3], comp[2]] == -2 ? ((:F, 4), comp) : ((:F, 4), reverse(comp))
        (:G, 2) => gcm[comp[2], comp[1]] == -3 ? ((:G, 2), comp) : ((:G, 2), reverse(comp))

        # Non-simply laced affine types
        (:A, 1, :aff) => begin
            # If we have the (-2, -2) Cartan matrix, then it is symmetric so nothing more to do.
            gcm[comp[1], comp[2]] == -2 && return ((:A, 1, :aff), comp)

            # If we have the (-1, -4) Cartan matrix, the affine vertex should be the one with the -1 in its row.
            return (gcm[comp[2], comp[1]] == -1) ? ((:BC, 1, :aff), comp) : ((:BC, 1, :aff), reverse(comp))
        end

        (:B, n, :aff) => begin
            # We know n ≥ 3. We need to figure out whether the bond between n-1 and n is pointing to the
            # right (type B~n) or left (type B@n).
            return (gcm[comp[n], comp[n-1]] == -2) ? ((:B, n, :aff), comp) : ((:B, n, :dualaff), comp)
        end

        (:C, n, :aff) => begin
            # We know n ≥ 2, and Coxeter graph:   (n+1) == 1 -- 2 -- ... -- (n-1) == n.

            # Both arrows pointing in is type C~n.
            gcm[comp[1], comp[n+1]] == -2 && gcm[comp[n-1], comp[n]] == -2 && return ((:C, n, :aff), comp)

            # Both arrows pointing out is type C@n.
            gcm[comp[1], comp[n+1]] == -1 && gcm[comp[n-1], comp[n]] == -1 && return ((:C, n, :dualaff), comp)

            # Both arrows pointing the same way is type BC~n. The two arrows both point away from the affine vertex.
            # (n+1) =>= 1 -- 2 ... -- (n-1) =>= n.
            gcm[comp[1], comp[n+1]] == -2 && return ((:BC, n, :aff), comp)

            # At this point we have the right diagram but in the wrong direction:
            # (n+1) =<= 1 -- 2 ... -- (n-1) =<= n.
            # The correct order is (n-1), (n-2), ..., 1, n+1, n.
            return ((:BC, n, :aff), [comp[n-1:-1:1]; comp[n+1]; comp[n]])
        end

        (:F, n, :aff) => begin
            # Arrow between 2 and 3 pointing right is usual affinisation, otherwise dual affinisation.
            return (gcm[comp[3], comp[2]] == -2) ? ((:F, 4, :aff), comp) : ((:F, 4, :dualaff), comp)
        end

        (:G, n, :aff) => begin
            # Arrow between 1 and 2 pointing right is usual affinisation, otherwise dual.
            return (gcm[comp[2], comp[1]] == -3) ? ((:G, 2, :aff), comp) : ((:G, 2, :dualaff), comp)
        end
    end

    return [convert_component(comp...) for comp in cox_components]
end

# Convert a type (:A, 4) to a string like "A4".
type_to_string(type::Tuple) = @match type begin
    (letter, rank, :aff) => "$letter~$rank"
    (letter, rank, :dualaff) => "$letter@$rank"
    (:I, 2, m) => "I2(m)"
    _ => join(map(string, type))
end

# Convert a list of types to a string like "A4 x H3"
type_to_string(types::Vector) = join([type_to_string(type) for (type, comp) in types], " x ")

# Check if a particular type like (:A, 5) is finite type.
is_finite_type(type::Tuple) = @match type begin
    (:Unknown, n) => false
    (:A, n) => true
    (:B, n) => true
    (:C, n) => true
    (:D, n) => true
    (:E, n) => 6 <= n <= 8
    (:F, 4) => true
    (:G, 2) => true
    (:H, n) => 2 <= n <= 4
    (:I, 2, m) => m >= 2
    _ => false
end

# Retrieve the degrees of a finite type Coxeter group.
degrees(type::Tuple) = @match type begin
    (:A, n) => Vector(2:n+1)
    (:B, n) => Vector(2:2:2*n)
    (:C, n) => Vector(2:2:2*n)
    (:D, n) => sort([Vector(2:2:2*n-2); n])
    (:E, 6) => [2, 5, 6, 8, 9, 12]
    (:E, 7) => [2, 6, 8, 10, 12, 14, 18]
    (:E, 8) => [2, 8, 12, 14, 18, 20, 24, 30]
    (:F, 4) => [2, 6, 8, 12]
    (:G, 2) => [2, 6]
    (:H, 3) => [2, 6, 10]
    (:H, 4) => [2, 12, 20, 30]
    (:I, 2, m) => [2, m]
    _ => error("$(type) is not finite type")
end

# Check if a list of types like [(:A, 5), (:B, 3, :aff)] is finite type.
is_finite_type(composite_type::Vector) = all(is_finite_type(type) for (type, comp) in composite_type)

# Check if a particular type like (:A, 5, :aff) is affine type.
is_affine_type(type::Tuple) = @match type begin
    (_, _, :aff) => true
    (_, _, :dualaff) => true
    _ => false
end

# Check if a list of types is affine type. This only occurs when there is a single type, and it is affine.
is_affine_type(composite_type::Vector) = length(composite_type) == 1 && is_affine_type(composite_type[1][1])


# A Coxeter system is a classified Coxeter matrix, along with some extra data.
struct CoxeterSystem
    coxeter_matrix::Matrix{Int}
    components::Vector{Tuple{Tuple, Vector{Int}}}
end

"""
    coxeter_system(mat)

Create a Coxeter system from a Coxeter matrix or Cartan matrix.
"""
function coxeter_system(mat::AbstractMatrix{T}) where {T <: Integer}
    if is_gcm(mat)
        mat = gcm_to_coxeter_matrix(mat)
    end
    is_coxeter_matrix(mat) || error("Given matrix was not a Coxeter or Cartan matrix.")
    components = classify_coxeter_matrix(mat)
    return CoxeterSystem(mat, components)
end

# Pretty-print some data about a Coxeter system.
function Base.show(io::IO, cox::CoxeterSystem)
    adjectives = [
        "rank $(rank(cox))",
        length(cox.components) == 1 ? "irreducible" : "reducible",
        is_finite_type(cox.components) ? "finite type" : is_affine_type(cox.components) ? "affine type" : "indefinite type",
    ]
    write(io, "Coxeter system ($(join(adjectives, ", "))) of type $(type_to_string(cox.components))")
end

"""
    coxeter_matrix(cox::CoxeterSystem)

Return the underlying Coxeter matrix of a Coxeter system.
"""
coxeter_matrix(cox::CoxeterSystem) = copy(cox.coxeter_matrix)

"""
    rank(cox::CoxeterSystem)

The rank of a Coxeter system is the number of simple generators.
"""
rank(cox::CoxeterSystem) = size(cox.coxeter_matrix)[1]

"""
    is_finite_type(cox::CoxeterSystem)

Return true if the Coxeter system is finite type, meaning that the associated Coxeter group is finite.
"""
is_finite_type(cox::CoxeterSystem) = is_finite_type(cox.components)

"""
    degrees(cox::CoxeterSystem)

Retrieve the degrees of a finite-type Coxeter system, or error if the system is not finite type.
"""
function degrees(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return sort([degree for (type, comp) in cox.components for degree in degrees(type)])
end

"""
    exponents(cox::CoxeterSystem)

The exponents of a finite-type Coxeter system (one less than the degrees), or error if the system is not finite type.
"""
function exponents(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return (degrees(cox) .- 1)
end

"""
    order(cox::CoxeterSystem)

Return the order of the group for a finite-type Coxeter system, or error if the system is not finite type.
"""
function order(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return reduce(*, degrees(cox); init=BigInt(1))
end

"""
    number_of_reflections(cox::CoxeterSystem)

Return the number of reflections (number of positive roots) in a finite-type Coxeter system, or error if the system is
not finite-type.
"""
function number_of_reflections(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return sum(degrees(cox)) - rank(cox)
end

"""
    coxeter_number(cox::CoxeterSystem)

Return the Coxeter number of the system, the order of the Coxeter element, or error if the system is not finite type.
"""
function coxeter_number(cox::CoxeterSystem)
    is_finite_type(cox) || error("Coxeter system has infinite order")
    return degrees(cox)[end]
end



# Get a Cartan matrix, in the Kac numbering convention.
cartan_matrix(type::Tuple) = @match type begin
    # Finite type GCMs.
    (:A, n), if n >= 0 end => [i == j ? 2 : abs(i - j) == 1 ? -1 : 0 for i in 1:n, j in 1:n]
    (:B, n), if n >= 2 end => begin
        gcm = cartan_matrix((:A, n))
        gcm[n, n-1] = -2
        return gcm
    end
    (:BC, 1) => [2;;]
    (:BC, n), if n >= 2 end => cartan_matrix((:B, n))
    (:C, n), if n >= 2 end => begin
        gcm = cartan_matrix((:A, n))
        gcm[n-1, n] = -2
        return gcm
    end
    (:D, n), if n >= 2 end => begin
        # Start with a block containing a path and two disconnected vertices, then connect them.
        gcm = cat(cartan_matrix((:A, n-2)), [2 0 ; 0 2]; dims=(1, 2))
        if n >= 3
            gcm[n, n-2] = gcm[n-2, n] = -1
            gcm[n, n-1] = gcm[n-1, n] = -1
        end
        return gcm
    end
    (:E, n), if 6 <= n <= 8 end => begin
        # Start with a path and a disconnected vertex, then attach the vertex at the right point.
        # Under the Kac convention, the vertex attaches to 3 in E6 and E7, and 5 in E8.
        mat = cat(cartan_matrix((:A, n-1)), [2;;]; dims=(1, 2))
        tri = (n == 6 || n == 7) ? 3 : 5
        mat[tri, n] = mat[n, tri] = -1
        return mat
    end
    (:F, 4) => [
        2 -1  0  0
       -1  2 -1  0
        0 -2  2 -1
        0  0 -1  2
    ]
    (:G, 2) => [
         2 -1
        -3  2
    ]

    # For affine type GCMs, build the root system and perform affinisation or dual affinisation.
    (label, n, :aff) => begin
        data = finite_gcm_data((label, n))
        data === nothing && error("The type $((label, n)) is not a finite type, and cannot be affinised.")
        data.is_irreducible || error("The finite type $((label, n)) is not irreducible, and cannot be affinised.")

        # To affinise, we want to add a root γ such that ⟨α_i^, γ⟩ = ⟨α_i^, -α~⟩, where α~ is the highest root in the
        # original system. Likewise, the new coroot γ^ should satisfy ⟨γ^, α_i⟩ = ⟨-(α~)^, α_i⟩, where (α~)^ is the
        # root dual to the highest root, often called the "highest short coroot". If the original Cartan matrix is A,
        # it then follows that the affinisation is
        # [ A            A*(-α~) ]
        # [ -(α~)^ * A   2       ]
        gcm = cartan_matrix((label, n))
        return [
            gcm                                             -gcm * reshape(data.highest_root, :, 1)
            -reshape(data.highest_root_dual, 1, :) * gcm     2
        ]
    end
    (label, n, :dualaff) => begin
        data = finite_gcm_data((label, n))
        data === nothing && error("The type $((label, n)) is not a finite type, and cannot be affinised.")
        data.is_irreducible || error("The finite type $((label, n)) is not irreducible, and cannot be affinised.")

        # Follow the comment above for :aff, however here the roles of roots and coroots are switched, and the GCM is
        # transposed (since we are really building an affine Weyl group as a subgroup of affine transformations in the
        # root space, and so the "roots" of the new system live in the affinised coroot space).
        gcm = transpose(cartan_matrix((label, n)))
        return [
            gcm                                               -gcm * reshape(data.highest_coroot, :, 1)
            -reshape(data.highest_coroot_dual, 1, :) * gcm     2
        ]
    end
    _ => error("No known Cartan matrix associated to type $(type)")
end

# A root or coroot in a root system
struct Root
    index::Int          # Index in the root system
    rt::Vector{Int}     # Coordinates in root basis
    wt::Vector{Int}     # Coordinates in weight basis
end

# A finite root system
struct RootSystem
    # Underlying Cartan matrix and divisibility vector (all 1's unless there are type BC components).
    gcm::Matrix{Int}
    divs::Vector{Int}

    # Type and the components giving each type.
    types::Vector{Tuple}
    comps::Vector{Vector{Int}}

    # Roots and coroots in parallel arrays, so the α ↦ α^ bijection is looking up the coroot with the same index.
    roots::Vector{Root}
    coroots::Vector{Root}

    # Look up a root by its coordinates in the root basis, or a coroot by its coordinates in the coroot basis.
    rtToIndex::Dict{Vector{Int}, Int}
    cortToIndex::Dict{Vector{Int}, Int}
end



"""
    root_system(gcm)

Construct a root system for the finite-type generalised Cartan matrix gcm.
"""
root_system(gcm::AbstractMatrix) = root_system(gcm, Set{Int}())

"""
    root_system(gcm, doubled)

Construct a root system for the finite-type generalised Cartan matrix gcm, with the specified simple roots doubled
to produce a non-reduced root system.
"""
function root_system(gcm::AbstractMatrix{T}, doubled::Set{Int}) where {T <: Integer}
    # Classify the types and components of the GCM, and check we are starting with a finite-type GCM.
    components = classify_gcm(gcm)
    is_finite_type(components) || error("Can only build the root system of a finite-type GCM.")
    types = [type for (type, comp) in components]
    comps = [comp for (type, comp) in components]

    # Ensure the doubled roots are in compatible places.
    rank = size(gcm)[1]
    if !all(1 <= i <= rank && all(gcm[i, j] % 2 == 0 for j in 1:rank) for i in doubled)
        error("Doubled roots may only be placed where the rows of the Cartan matrix are divisible by two.")
    end

    # Create the divided Cartan matrix dgcm = D^-1 A, and a vector of 1's and 2's called divs.
    divs = [i ∈ doubled ? 2 : 1 for i in 1:rank]
    # print("Divs: $(divs)\n")
    dgcm = copy(gcm)
    for i in doubled
        dgcm[i, :] .÷= 2
    end

    # Create the initial (coroot, root) pairs for the simple roots.
    roots = [Root(i, [i == j ? 1 : 0 for j in 1:rank], dgcm[:, i]) for i in 1:rank]
    coroots = [Root(i, [i == j && i ∈ doubled ? 2 : i == j ? 1 : 0 for j in 1:rank], gcm[i, :]) for i in 1:rank]

    # Add the (coroot, root) pairs for the doubled roots.
    for i in doubled
        push!(roots, Root(length(roots) + 1, [i == j ? 2 : 0 for j in 1:rank], 2 * gcm[:, i]))
        push!(coroots, Root(length(coroots) + 1, [i == j ? 1 : 0 for j in 1:rank], gcm[i, :] .÷ 2))
    end

    # Add the coroots and roots we've seen so far to lookup tables.
    rtToIndex = Dict(root.rt => root.index for root in roots)
    cortToIndex = Dict(coroot.rt => coroot.index for coroot in coroots)

    # We'll need a function which takes a (root, coroot) pair and applies the simple reflection s.
    function reflect_pair(root::Root, coroot::Root, s::Int)
        # s(root) = root - ⟨α_s^, root⟩ α_s. The pairing ⟨α_s^, root⟩ is d_s wt[s]. In the root basis this subtracts
        # a single coordinate, and in the weight basis this subtracts a vector.
        pairing = divs[s] * root.wt[s]
        # print("Pairing of rt=$(root.rt), wt=$(root.wt) with α_$(s)^ is $(pairing)\n")
        new_root = Root(root.index + 1, copy(root.rt), copy(root.wt))
        new_root.rt[s] -= pairing
        new_root.wt .-= pairing * dgcm[:, s]

        # s(coroot) = coroot - ⟨coroot, α_s⟩ α_s^. The copairing ⟨coroot, α_s⟩ is wt[s].
        copairing = coroot.wt[s]
        new_coroot = Root(coroot.index + 1, copy(coroot.rt), copy(coroot.wt))
        new_coroot.rt[s] -= copairing * divs[s]
        new_coroot.wt .-= copairing * gcm[s, :]

        return (new_root, new_coroot)
    end

    # Convenience function for adding a pair of roots to the data structures, or ignoring them if they exist already.
    function add_pair(root::Root, coroot::Root)
        if root.rt ∉ keys(rtToIndex)
            push!(roots, root)
            push!(coroots, coroot)
            rtToIndex[root.rt] = root.index
            cortToIndex[coroot.rt] = coroot.index
        end
    end

    # Now play the root reflection game: keep applying reflections, making a breadth-first search of the root
    # reflection graph, avoiding the negative roots. This gives us all the positive roots.
    pos = 1
    while pos < length(roots)
        (root, coroot) = (roots[pos], coroots[pos])
        pos += 1

        for s in 1:rank
            (sroot, scoroot) = reflect_pair(root, coroot, s)
            # print("Reflection of rt=$(root.rt) wt=$(root.wt) in $(s) is rt=$(sroot.rt), wt=$(sroot.wt)\n")

            # Ignore negative roots
            if all(x -> x >= 0, sroot.rt)
                add_pair(sroot, scoroot)
            # else
                # print("Discarding $(root.rt)\n")
            end
        end
    end

    # Add all negative roots.
    nPositive = length(roots)
    for i in 1:nPositive
        add_pair(Root(nPositive + i, -roots[i].rt, -roots[i].wt), Root(nPositive + i, -coroots[i].rt, -coroots[i].wt))
    end

    return RootSystem(gcm, divs, types, comps, roots, coroots, rtToIndex, cortToIndex)
end

function Base.show(io::IO, rs::RootSystem)
    write(io, "Root system of type $(rs.types), rank $(rank(rs)) with $(length(rs.roots)) roots.")
end

"""
    rank(rs::RootSystem)

Return the rank of the root system, the number of simple roots.
"""
rank(rs::RootSystem) = size(rs.gcm)[1]

"""
    highest_root(rs::RootSystem; basis=:root)

Return the highest root, in the root basis. Pass basis=:weight for the weight basis.
"""
function highest_root(rs::RootSystem; basis=:root)
    (basis == :root || basis == :weight) || error("The basis must be :root or :weight")

end

"""
    positive_roots(rs::RootSystem, basis=:root)

Return the positive roots, in the root basis. Pass basis:=weight for the weight basis.
"""
function positive_roots(rs::RootSystem; basis=:root)
    (basis == :root || basis == :weight) || error("The basis must be :root or :weight")
    return hcat([(basis == :root) ? x.rt : x.wt for x in rs.roots]...)
end

"""
    positive_coroots(rs::RootSystem, basis=:root)

Return the positive coroots, in the coroot basis. Pass basis:=weight for the coweight basis.
"""
function positive_coroots(rs::RootSystem; basis=:root)
    (basis == :root || basis == :weight) || error("The basis must be :root or :weight")
    return hcat([(basis == :root) ? x.rt : x.wt for x in rs.coroots[1:length(rs.coroots)÷2]]...)
end



# This type should not be used outside of this file. It is used to construct the affinisations and dual affinisations
# of finite-type irreducible root systems.
#
# This is giving the wrong thing when affinising for type BC because the pairing between the root and coroot bases
# is no longer the Cartan matrix A, instead it's D^-1 A.
# Perhaps it would be better after all to build the root system, and then just take the fundamental weight vectors
# for the highest root/coroot etc.
struct FiniteGCMData
    is_irreducible::Bool                # Is this type irreducible (so it has an affinisation).
    highest_root::Vector{Int}           # Highest root, in the root basis.
    highest_root_dual::Vector{Int}      # (Highest root)^, in the coroot basis. Also called the highest short coroot.
    highest_coroot::Vector{Int}         # Highest coroot, in the coroot basis.
    highest_coroot_dual::Vector{Int}    # (Highest coroot)^, in the root basis. Also called the highest short root.
end

# Shortcut method for simply-laced types where highest roots are all the same.
FiniteGCMData(is_irreducible::Bool, root::Vector{Int}) = FiniteGCMData(is_irreducible, root, root, root, root)

# Given a finite-type label, return a FiniteGCMData giving the highest roots,
# or nothing if the system is not defined.
finite_gcm_data(type::Tuple) = @match type begin
    (:A, n) => FiniteGCMData(n >= 1, [1 for i in 1:n])
    (:B, n) => FiniteGCMData(
        n >= 2,
        [i == 1 ? 1 : 2 for i in 1:n],
        [i == n ? 1 : 2 for i in 1:n],
        [1 for i in 1:n],
        [(i == 1 || i == n) ? 1 : 2 for i in 1:n],
    )
    (:BC, n) => FiniteGCMData(
        n >= 1,
        [2 for i in 1:n],
        [1 for i in 1:n],
        [2 for i in 1:n],
        [1 for i in 1:n],
    )
    (:C, n) => FiniteGCMData(
        n >= 2,
        [i == n ? 1 : 2 for i in 1:n],
        [i == 1 ? 1 : 2 for i in 1:n],
        [(i == 1 || i == n) ? 1 : 2 for i in 1:n],
        [1 for i in 1:n],
    )
    (:D, n) => FiniteGCMData(n >= 3, [(i == 1 || i == n-1 || i == n) ? 1 : 2 for i in 1:n])
    (:E, 6) => FiniteGCMData(true, [1, 2, 3, 2, 1, 2])
    (:E, 7) => FiniteGCMData(true, [2, 3, 4, 3, 2, 1, 2])
    (:E, 8) => FiniteGCMData(true, [2, 3, 4, 5, 6, 4, 2, 3])
    (:F, 4) => FiniteGCMData(
        true,
        [2, 3, 4, 2],
        [2, 4, 3, 2],
        [1, 2, 3, 2],
        [2, 3, 2, 1],
    )
    (:G, 2) => FiniteGCMData(true, [2, 3], [3, 2], [1, 2], [2, 1])
    _ => nothing
end


# (0, 2)=[-2,-2] (1, 2)        (2, 2) = s1(2a2)
#        (0, 1)          (1, 1) = s1(a2)
#                  0           (1, 0) = [2, 0]
#

# 2 ([2, 0], [-2, -2])
#  -----------------    = -8/4 = -2
# ([2, 0], [2, 0])
