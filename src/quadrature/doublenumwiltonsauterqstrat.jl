function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumWiltonSauterQStrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule_far,qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far,qs.inner_rule_near))
     
    leg = (
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_vert,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_edge,0,1)),
      convert.(NTuple{2,T},_legendre(qs.sauter_schwab_common_face,0,1)),)

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end


function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace,  i, τ, j, σ, qd,
    qs::DoubleNumWiltonSauterQStrat)

    T = eltype(eltype(verticeslist(τ)))
    hits = 0
    dtol = 1.0e3 * eps(T)
    dmin2 = floatmax(T)
    for t in verticeslist(τ)
        for s in verticeslist(σ)
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end

    hits == 3 && return SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    hits == 2 && return SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    hits == 1 && return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))
    if max(dmin2*k2, dmin2/16h2) < xtol2
        return WiltonSERule(
            qd.tpoints[2,i],
            DoubleQuadRule(
                qd.tpoints[2,i],
                qd.bpoints[2,j],),)
    end

    return DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j],)
end