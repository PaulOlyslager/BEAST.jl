struct Identity <: LocalOperator
end

kernelvals(biop::Identity, x) = nothing
integrand(op::Identity, kernel, x,y, g, f) = _krondot(getvalue(g), getvalue(f))
scalartype(op::Identity) = Union{}

struct NCross <: LocalOperator
end

kernelvals(op::NCross, mp) = nothing
integrand(op::NCross, kernel, x,y, g, f) = _krondot(getvalue(g), Ref(normal(x)) .× getvalue(f))
scalartype(op::NCross) = Union{}

function _alloc_workspace(qd, g, f, tels, bels)
    q = qd[1]
    τ = tels[1]
    σ = bels[1]
    w, p1, p2 = q[1], neighborhood(τ,q[2]), neighborhood(σ,q[2])
    a = (w, (p1,p2), g(p1), f(p2))
    A = Vector{typeof(a)}(undef,length(qd))
end

const LinearRefSpaceTriangle = Union{RTRefSpace, NDRefSpace, BDMRefSpace, NCrossBDMRefSpace}
defaultquadstrat(::LocalOperator, ::LinearRefSpaceTriangle, ::LinearRefSpaceTriangle) = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::Union{LinearRefSpaceTriangle,LagrangeRefSpace{T,Deg,3}}, f::Union{LagrangeRefSpace,LinearRefSpaceTriangle}, tels, bels,
        qs::SingleNumQStrat) where{T,Deg}

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end

defaultquadstrat(::LocalOperator, ::subReferenceSpace, ::subReferenceSpace) = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::subReferenceSpace, f::subReferenceSpace, tels, bels,
        qs::SingleNumQStrat)

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

const LinearRefSpaceTetr = Union{NDLCCRefSpace, NDLCDRefSpace, BDM3DRefSpace}
defaultquadstrat(::LocalOperator, ::LinearRefSpaceTetr, ::LinearRefSpaceTetr) = SingleNumQStrat(3)
function quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr, tels, bels, qs::SingleNumQStrat)

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = _alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,2}, ::LagrangeRefSpace{T,D2,2}) where {T,D1,D2} = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,2} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)
    U = typeof(tels[1].volume)
    u, w = legendre(qs.quad_rule, 0.0, 1.0)
    qd = [(w[i],u[i]) for i in eachindex(w)]
    A = _alloc_workspace(Tuple{U,U}.(qd), g, f, tels, bels)
    return qd, A
end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,3}, ::LagrangeRefSpace{T,D2,3}) where {T,D1,D2} = SingleNumQStrat(6)
# function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,3} where {T,Deg},
#     f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
#     A = _alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end

defaultquadstrat(::LocalOperator, ::LagrangeRefSpace{T,D1,4}, ::LagrangeRefSpace{T,D2,4}) where {T,D1,D2} = SingleNumQStrat(6)
function quaddata(op::LocalOperator, g::LagrangeRefSpace{T,Deg,4} where {T,Deg},
    f::LagrangeRefSpace, tels::Vector, bels::Vector, qs::SingleNumQStrat)

     o, x, y, z = CompScienceMeshes.euclidianbasis(3)
     reftet = simplex(x,y,z,o)
     qps = quadpoints(reftet, qs.quad_rule)
     qd = [(w, parametric(p)) for (p,w) in qps]
     A = _alloc_workspace(qd, g, f, tels, bels)
     return qd, A
end


function quadrule(op::LocalOperator, ψ::RefSpace, ϕ::RefSpace, τ,σ, (qd,A), qs::SingleNumQStrat)
    for i in eachindex(qd)
        q = qd[i]
        w, p1, p2 = q[1], neighborhood(τ,q[2]), neighborhood(σ,q[2])
        A[i] = (w, (p1,p2), ψ(p1), ϕ(p2))

    end
    return A
end

# function quadrule(op::HHHLocalOperator, ψ::RefSpace, ϕ::RefSpace, τ, (qd,A),nt,nb, qs::SingleNumQStrat)
#     for i in eachindex(qd)
#         q = qd[i]
#         w, p = q[1], extendedneighborhood(τ,q[2],nt,nb)
#         A[i] = (w, p, ψ(p), ϕ(p))
#     end
#     return A
# end


