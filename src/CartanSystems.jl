export cartan_matrix, root_system, positive_roots, positive_coroots

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


# Return a (2 x rank) matrix of (x, y) coordinates for each vertex, to lay out the graph.
function layout_coxeter_nodes(mat::AbstractMatrix{T}) where {T <: Integer}
    coxmat = is_gcm(mat) ? gcm_to_coxeter_matrix(mat) : mat
    is_coxeter_matrix(coxmat) || error("Matrix was not a GCM or Coxeter matrix")

    # Adjacency matrix of underlying undirected graph.
    rank = size(mat)[1]
    adj = Dict(i => [j for j in 1:rank if i != j && coxmat[i, j] != 2] for i in 1:rank)
    deg = [length(adj[i]) for i in 1:rank]

    # Perform a breadth-first search from a point, returning a triple (vertices, dist, pred) of the vertices visited
    # in BFS order, their distances from the source, and a dictionary mapping each vertex to its predecessor on a
    # shortest path to the source.
    function bfs(adj, start)
        dist = Dict(start => 0)
        pred = Dict(start => start)
        order = [start]
        pos = 1
        while pos <= length(order)
            s = order[pos]
            pos += 1
            for t in adj[s]
                if !haskey(dist, t)
                    dist[t] = dist[s] + 1
                    pred[t] = s
                    push!(order, t)
                end
            end
        end
        return (order, dist, pred)
    end

    # Given a vertex in a cycle, return a list of vertices in that cycle, starting from start, in cycle order.
    function getcycle(adj, start)
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

    # Lay out a single connected component.
    function layout_component(comp::Vector{Int})
        # A connected component is a cycle iff every vertex has degree 2.
        if all([deg[s] == 2 for s in comp])
            order = getcycle(adj, minimum(comp))

    end

end
