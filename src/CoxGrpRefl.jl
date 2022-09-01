export CoxGrpRefl, CoxEltRefl, coxeter_group_reflection

using LinearAlgebra: I

mutable struct CoxGrpRefl <: CoxGrp
    cartan_mat::Matrix{Int64}
end

struct CoxEltRefl <: CoxElt
    group::CoxGrpRefl
    fwd::Matrix{Int64}
    inv::Matrix{Int64}
end

"""
    coxeter_group_reflection(gcm::Matrix) -> CoxGrpRefl

Given a GCM, create a Coxeter group whose elements are represented by matrices in the reflection representation.
"""
function coxeter_group_reflection(gcm::Matrix{T}) where T <: Integer
    is_gcm(gcm) || error("The argument must be a GCM.")
    grp = CoxGrpRefl(convert(Matrix{Int64}, gcm))
    return grp, generators(grp)
end

rank(grp::CoxGrpRefl) = size(grp.cartan_mat)[1]

function Base.one(grp::CoxGrpRefl)
    rank = size(grp.cartan_mat)[1]
    id = Matrix{Int64}(I, rank, rank)
    return CoxEltRefl(grp, id, id)
end

function generators(grp::CoxGrpRefl)
    rank = size(grp.cartan_mat)[1]
    gens = [Matrix{Int64}(I, rank, rank) for _ in 1:rank]
    for i in 1:rank
        gens[i][:, i] .-= grp.cartan_mat[:, i]
    end
    return [CoxEltRefl(grp, gen, gen) for gen in gens]
end

Base.parent(w::CoxEltRefl) = w.group
