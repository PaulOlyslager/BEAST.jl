#### implementation of the user defined quadstrat. wrappers can be written to map to the existing quadstrats.
#### the idea is to overcome the cumbersome code repeating hussle of copieing all existing rules if a new ruels is defiend to create a new quadstrat.
#### dict makes it slow, reserve the memory beforehand. (maybe move it to quaddata)
abstract type RefQuadStrat <: AbstractQuadStrat end
abstract type SelectionCriteria end
#TODO check if nested momentitegrals in quadstrat is faster than dispatch afterwards
"""
    Creates the dedicated Rule for the given information. cellinfo etc is given in kwargs.
"""
function (t::RefQuadStrat)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,
    quad_data , zlocal, test_space, tptr, trial_space, bptr)
    display(t)
    error("quadrule not implemented")
end
# function quadrule(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,
#     quad_data, qs::RefQuadStrat)


#     qs(operator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,nothing,quad_data)
# end
function quadrule_momintegrals!(zlocal, biop,
            test_space,  tptr, tcell,
            trial_space, bptr, bcell, test_shapes, trial_shapes, p, q, qd, qs)

    return qs(biop,
    test_shapes, trial_shapes,
    p, tcell, q, bcell,nothing,qd, zlocal, test_space, tptr, trial_space, bptr)
    
end
function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qs::RefQuadStrat)
    error("quaddata not implemented")
end
function quaddatathreadlocal(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qd, qs::RefQuadStrat)
    qd
end
"""
    returns a boolean and a list of kwargs.
"""
function (::SelectionCriteria)(operator,
            local_test_basis, local_trial_basis,
            test_id, test_element, trial_id, trial_element,cache)
    error("selection criteria not implemented")
end

struct CombinedRefQuadStrat{U, V} <: RefQuadStrat
    rules::U
    criteria::V
end
function (cr::CombinedRefQuadStrat)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache, qd )
    _setzero!(qd.cache)
    for (rule, crit, qdit) in zip(cr.rules, cr.criteria,qd.qd)
        val = crit(operator,
            local_test_basis, local_trial_basis,
            test_id, test_element, trial_id, trial_element,qd.cache) 
        if val
            return rule(operator,
            local_test_basis, local_trial_basis,
            test_id, test_element, trial_id, trial_element,qd.cache,qdit)
        end
    end
    error("no rule selected")
end
# function (cr::CombinedRefQuadStrat)(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache, qd,  
#     zlocal, test_space, tptr, trial_space, bptr )
#     _setzero!(qd.cache)
#     for (rule, crit, qdit) in zip(cr.rules, cr.criteria,qd.qd)
#         val = crit(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache) 
#         if val
#             rule(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache,qdit,
#             zlocal, test_space, tptr, trial_space, bptr)
#             return nothing
#         end
#     end
#     error("no rule selected")
# end

# function (cr::CombinedRefQuadStrat{<:Tuple{U1,U2,U3,U4},V})(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache, qd,  
#     zlocal, test_space, tptr, trial_space, bptr ) where {U1, U2, U3, U4, V}
#     # return cr, operator, local_test_basis, local_trial_basis, test_id, test_element, trial_id, trial_element,cache, qd, zlocal, test_space, tptr, trial_space, bptr
#     _setzero!(qd.cache)
#     r1, r2, r3, r4 = cr.rules[1], cr.rules[2], cr.rules[3], cr.rules[4]
#     c1, c2, c3, c4 = cr.criteria[1], cr.criteria[2], cr.criteria[3], cr.criteria[4]
#     qd1, qd2, qd3, qd4 = qd.qd[1], qd.qd[2], qd.qd[3], qd.qd[4]
#     val = c1(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache)
#     if val
#         r1(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache,qd1,
#             zlocal, test_space, tptr, trial_space, bptr)
#         return nothing
#     end
#         val = c2(operator,
#                 local_test_basis, local_trial_basis,
#                 test_id, test_element, trial_id, trial_element,qd.cache)
#     if val
#         r2(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache,qd2,
#             zlocal, test_space, tptr, trial_space, bptr)
#         return nothing
#     end
#         val = c3(operator,
#                 local_test_basis, local_trial_basis,
#                 test_id, test_element, trial_id, trial_element,qd.cache)
#     if val
#         r3(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache,qd3,
#             zlocal, test_space, tptr, trial_space, bptr)
#         return nothing
#     end
#         val = c4(operator,
#                 local_test_basis, local_trial_basis,
#                 test_id, test_element, trial_id, trial_element,qd.cache)
#     if val
#         r4(operator,
#             local_test_basis, local_trial_basis,
#             test_id, test_element, trial_id, trial_element,qd.cache,qd4,
#             zlocal, test_space, tptr, trial_space, bptr)
#         return nothing
#     end
#      error("no rule selected")
# end
@generated function (cr::CombinedRefQuadStrat{R,C})(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element, cache, qd,
    zlocal, test_space, tptr, trial_space, bptr) where {R<:Tuple, C}

    N = length(R.parameters)

    branches = []

    for i in 1:N
        push!(branches, quote
            if cr.criteria[$i](operator,
                local_test_basis, local_trial_basis,
                test_id, test_element, trial_id, trial_element, qd.cache)

                cr.rules[$i](operator,
                    local_test_basis, local_trial_basis,
                    test_id, test_element, trial_id, trial_element,
                    qd.cache, qd.qd[$i],
                    zlocal, test_space, tptr, trial_space, bptr)

                return nothing
            end
        end)
    end

    return quote
        _setzero!(qd.cache)
        $(Expr(:block, branches...))
        error("no rule selected")
    end
end
# @generated function (cr::CombinedRefQuadStrat{U,V})(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache, qd) where {U,V}
#     ex = quote
#         _setzero!(qd.cache)
#     end
#         for (rule, crit, qdit) in zip(cr.rules, cr.criteria,qd.qd)
#             val = crit(operator,
#                 local_test_basis, local_trial_basis,
#                 test_id, test_element, trial_id, trial_element,qd.cache) 
#             if val
#                 return rule(operator,
#                 local_test_basis, local_trial_basis,
#                 test_id, test_element, trial_id, trial_element,qd.cache,qdit)
#             end
#         end
#         error("no rule selected")
#     end
# end

function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, cr::CombinedRefQuadStrat) 
    Tuple(quaddata(operator, local_test_basis, local_trial_basis, test_elements, trial_elements, rule) for rule in cr.rules)
end
# keys = (:operator, :local_test_basis, :local_trial_basis, :test_elements, :trial_elements, :dist)
# function quaddatathreadlocal(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_elements, trial_elements, qd, cr::CombinedRefQuadStrat)
#     (qd = Tuple(quaddatathreadlocal(operator, local_test_basis, local_trial_basis, test_elements, trial_elements, qd[i], cr.rules[i]) for i in 1:length(cr.rules)),
#     cache = (dist = Ref{Union{Float64,Nothing}}(nothing),
#              hits = Ref{Union{Int64,Nothing}}(nothing),
#              idx_t = Ref{Union{Vector{Int64},Nothing}}(nothing),
#              idx_s = Ref{Union{Vector{Int64},Nothing}}(nothing)))
# end
@generated function quaddatathreadlocal(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qd, cr::CombinedRefQuadStrat{R,U}) where {R<:Tuple,U}

    N = length(R.parameters)

    # build tuple elements
    elems = []
    for i in 1:N
        push!(elems, :(quaddatathreadlocal(
            operator, local_test_basis, local_trial_basis,
            test_elements, trial_elements,
            qd[$i], cr.rules[$i]
        )))
    end

    return quote
        (
            qd = tuple($(elems...)),
            cache = (
                dist = Ref{Union{Float64,Nothing}}(nothing),
                hits = Ref{Union{Int64,Nothing}}(nothing),
                idx_t = Ref{Union{Vector{Int64},Nothing}}(nothing),
                idx_s = Ref{Union{Vector{Int64},Nothing}}(nothing)
            )
        )
    end
end
function _setzero!(cache) 
    cache.dist.x = nothing
    cache.hits.x = nothing
    cache.idx_t.x = nothing
    cache.idx_s.x = nothing
end

###### list of Selection Criteria
struct CRTrue <: SelectionCriteria end
function (cr::CRTrue)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    return true
end
struct CRWithinDist{R} <: SelectionCriteria
    distmin::R
    distmax::R
end
function (cr::CRWithinDist)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    if cache.dist[1] !== nothing
        dist = cache.dist.x
    else
        dist = _distance(test_element, trial_element)
        cache.dist.x = dist
    end
    return cr.distmin <= dist <= cr.distmax
end
_distance(testelement, trialelement) = minimum([norm(t-s) for t in testelement.vertices for s in trialelement.vertices])
struct CRCommonVertex <: SelectionCriteria end
struct CRCommonEdge <: SelectionCriteria end
struct CRCommonFace <: SelectionCriteria end
struct CRCommonVolume <: SelectionCriteria end
function (cr::CRCommonVertex)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    if cache.hits.x === nothing 
        hits = _hitinfo(test_element, trial_element)
            cache.hits.x = hits
            # cache.idx_t.x = idx_t
            # cache.idx_s.x = idx_s
    else
        hits = cache.hits.x
    end 
    return hits == 1
end
function (cr::CRCommonEdge)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    if cache.hits.x === nothing 
        hits = _hitinfo(test_element, trial_element)
            cache.hits.x = hits
            # cache.idx_t.x = idx_t
            # cache.idx_s.x = idx_s
    else
        hits = cache.hits.x
    end 
    return hits == 2
end
function (cr::CRCommonFace)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    if cache.hits.x === nothing 
        hits = _hitinfo(test_element, trial_element)
            cache.hits.x = hits
            # cache.idx_t.x = idx_t
            # cache.idx_s.x = idx_s
    else
        hits = cache.hits.x
    end 
    return hits == 3
end
function (cr::CRCommonVolume)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache)
    if cache.hits.x === nothing 
        hits = _hitinfo(test_element, trial_element)
            cache.hits.x = hits
            # cache.idx_t.x = idx_t
            # cache.idx_s.x = idx_s
    else
        hits = cache.hits.x
    end 
    return hits == 4
end
# function (cr::CRCommonEdge)(;test_element, trial_element, hitinfo=_hitinfo(test_element,trial_element), kwargs...)
#     return (hitinfo.hits == 2, (hitinfo=hitinfo, test_element, trial_element,kwargs...))
# end
# function (cr::CRCommonFace)(;test_element, trial_element, hitinfo=_hitinfo(test_element,trial_element), kwargs...)
#     return (hitinfo.hits == 3, (hitinfo=hitinfo, test_element, trial_element,kwargs...))
# end
# function (cr::CRCommonVolume)(;test_element, trial_element, hitinfo=_hitinfo(test_element,trial_element), kwargs...)
#     return (hitinfo.hits == 4, (hitinfo=hitinfo, test_element, trial_element,kwargs...))
# end
function _hitinfo(testelement, trialelement)
    dtol = 1.0e3 * eps(eltype(eltype(testelement.vertices)))

    hits = 0
    # idx_t = Int64[]
    # idx_s = Int64[]
    # sizehint!(idx_t,4)
    # sizehint!(idx_s,4)
    # dmin2 = floatmax(eltype(eltype(testelement.vertices)))
    D = dimension(testelement)+dimension(trialelement)
    for (i,t) in enumerate(testelement.vertices)
        for (j,s) in enumerate(trialelement.vertices)
            # d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            # dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                # push!(idx_t,i)
                # push!(idx_s,j)
                hits +=1
                break
            end
        end
    end
    return hits#idx_t=idx_t, idx_s=idx_s)
end

###### list of quadrature rules
struct RefDoubleNumRule{R} <: RefQuadStrat
    outer_rule::R 
    inner_rule::R
end
# function (rule::RefDoubleNumRule)(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache, qd)
#     DoubleQuadRule(qd.tpoints[1, test_id], qd.bpoints[1, trial_id])
    
# end
function (rule::RefDoubleNumRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache, qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = DoubleQuadRule(qd.tpoints[1, test_id], qd.bpoints[1, trial_id])
    return momintegrals!(zlocal, operator,
    test_space,  tptr, test_element,
    trial_space, bptr, trial_element, qr)
end

struct RefSIMDDoubleNumRule{R} <: RefQuadStrat
    outer_rule::R 
    inner_rule::R
end
# function (rule::RefSIMDDoubleNumRule)(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache, qd)
#     SIMDDoubleQuadRule(qd.tpoints[1, test_id], qd.bpoints[1, trial_id], qd.cache)
# end
function (rule::RefSIMDDoubleNumRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache, qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = SIMDDoubleQuadRule(qd.tpoints[1, test_id], qd.bpoints[1, trial_id], qd.cache)
    return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, rule::RefSIMDDoubleNumRule) 
    # simdoperator = SIMDOperator(operator)
    derivative_needed = requires_derivative(operator)
    (tpoints = quadpoints_simd(local_test_basis, test_elements, (rule.outer_rule,),derivative_needed), 
     bpoints = quadpoints_simd(local_trial_basis, trial_elements, (rule.inner_rule,),derivative_needed))
end
function quaddatathreadlocal(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, qd, rule::RefSIMDDoubleNumRule) 
    g_cache = generate_cacheSIMDDoubleNum(operator, qd.tpoints, qd.bpoints)
    return (qd..., cache = g_cache)
end
struct RefWiltonSERule{R} <: RefQuadStrat
    outer_rule::R 
    inner_rule::R
end
# function (rule::RefWiltonSERule)(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache,qd;)
#     return WiltonSERule(
#         qd.tpoints[1, test_id],
#         DoubleQuadRule(
#             qd.tpoints[1, test_id],
#             qd.bpoints[1, trial_id],),)
# end
function (rule::RefWiltonSERule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = WiltonSERule(
        qd.tpoints[1, test_id],
        DoubleQuadRule(
            qd.tpoints[1, test_id],
            qd.bpoints[1, trial_id],),)
    return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, rule::Union{RefDoubleNumRule, RefWiltonSERule}) 
    (tpoints = quadpoints(local_test_basis, test_elements, (rule.outer_rule,)), 
     bpoints = quadpoints(local_trial_basis, trial_elements, (rule.inner_rule,)))
end
abstract type RefSauterschwab4DRule <: RefQuadStrat end
abstract type RefSauterschwab5DRule <: RefQuadStrat end
abstract type RefSauterschwab6DRule <: RefQuadStrat end
struct RefSauter4DCommonVertexRule{R} <: RefSauterschwab4DRule
    sauter_schwab_common_vert::R
end
struct RefSauter5DCommonVertexRule{R} <: RefSauterschwab5DRule
    sauter_schwab_common_vert::R
end
struct RefSauter6DCommonVertexRule{R} <: RefSauterschwab6DRule
    sauter_schwab_common_vert::R
end
struct RefSauter4DCommonEdgeRule{R} <: RefSauterschwab4DRule
    sauter_schwab_common_edge::R
end
struct RefSauter5DCommonEdgeRule{R} <: RefSauterschwab5DRule
    sauter_schwab_common_edge::R
end
struct RefSauter6DCommonEdgeRule{R} <: RefSauterschwab6DRule
    sauter_schwab_common_edge::R
end
struct RefSauter4DCommonFaceRule{R} <: RefSauterschwab4DRule
    sauter_schwab_common_face::R
end
struct RefSauter5DCommonFaceRule{R} <: RefSauterschwab5DRule
    sauter_schwab_common_face::R
end
struct RefSauter6DCommonFaceRule{R} <: RefSauterschwab6DRule
    sauter_schwab_common_face::R
end
struct RefSauter6DCommonVolumeRule{R} <: RefSauterschwab6DRule
    sauter_schwab_common_tetr::R
end
getvalue(rule::RefSauter4DCommonVertexRule) = rule.sauter_schwab_common_vert
getvalue(rule::RefSauter5DCommonVertexRule) = rule.sauter_schwab_common_vert
getvalue(rule::RefSauter6DCommonVertexRule) = rule.sauter_schwab_common_vert
getvalue(rule::RefSauter4DCommonEdgeRule) = rule.sauter_schwab_common_edge
getvalue(rule::RefSauter5DCommonEdgeRule) = rule.sauter_schwab_common_edge
getvalue(rule::RefSauter6DCommonEdgeRule) = rule.sauter_schwab_common_edge
getvalue(rule::RefSauter4DCommonFaceRule) = rule.sauter_schwab_common_face
getvalue(rule::RefSauter5DCommonFaceRule) = rule.sauter_schwab_common_face
getvalue(rule::RefSauter6DCommonFaceRule) = rule.sauter_schwab_common_face
getvalue(rule::RefSauter6DCommonVolumeRule) = rule.sauter_schwab_common_tetr
# function (rule::RefSauter4DCommonVertexRule)(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,cache,qd)
#     return SauterSchwabQuadrature.CommonVertex(qd.gausslegendre)
# end
function (rule::RefSauter4DCommonVertexRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = SauterSchwabQuadrature.CommonVertex(qd.gausslegendre)
    return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter5DCommonVertexRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
     qr = SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[3],qd.sing_qp[2]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
    return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[3],qd.sing_qp[2]))
end
function (rule::RefSauter6DCommonVertexRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwab3D.CommonVertex6D_S(SauterSchwab3D.Singularity6DPoint(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[3],qd.sing_qp[2]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter4DCommonEdgeRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
     zlocal, test_space, tptr, trial_space, bptr)
     qr = SauterSchwabQuadrature.CommonEdge(qd.gausslegendre)
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end

function (rule::RefSauter5DCommonEdgeRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter6DCommonEdgeRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = SauterSchwab3D.CommonEdge6D_S(SauterSchwab3D.Singularity6DEdge(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3],qd.sing_qp[4]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter4DCommonFaceRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwabQuadrature.CommonFace(qd.gausslegendre)
        return momintegrals!(zlocal, operator,
            test_space,  tptr, test_element,
            trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter4DCommonFaceRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
    qr = SauterSchwabQuadrature.CommonFace(qd.gausslegendre)
    return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter5DCommonFaceRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter6DCommonFaceRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwab3D.CommonFace6D_S(SauterSchwab3D.Singularity6DFace(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3],qd.sing_qp[4]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end
function (rule::RefSauter6DCommonVolumeRule)(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_id, test_element, trial_id, trial_element,cache,qd,
    zlocal, test_space, tptr, trial_space, bptr)
        qr = SauterSchwab3D.CommonVolume6D_S(SauterSchwab3D.Singularity6DVolume(cache.idx_t.x,cache.idx_s.x),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[4]))
     return momintegrals!(zlocal, operator,
        test_space,  tptr, test_element,
        trial_space, bptr, trial_element, qr)
end

function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, rule::Union{RefSauterschwab5DRule, RefSauterschwab6DRule})
    val = getvalue(rule)
    sing_qp = (SauterSchwab3D._legendre(val,0,1), 
               SauterSchwab3D._shunnham2D(val),
               SauterSchwab3D._shunnham3D(val),
               SauterSchwab3D._shunnham4D(val),)
    return (sing_qp = sing_qp,)
end
function quaddata(operator::IntegralOperator,
    local_test_basis, local_trial_basis,
    test_elements, trial_elements, rule::RefSauterschwab4DRule)
    T = coordtype(test_elements[1])
    leg = convert.(NTuple{2,T},_legendre(getvalue(rule),0,1))
    return (gausslegendre = leg,)
end
function assemblechunk_body!(biop,
        test_space, test_elements, test_assembly_data, test_cell_ptrs,
        trial_space, trial_elements, trial_assembly_data, trial_cell_ptrs,
        qdl, zlocal, store, quadstrat::T) where {T<:RefQuadStrat}

    test_shapes = refspace(test_space)
    trial_shapes = refspace(trial_space)

    # qdl = quaddatathreadlocal(biop, test_shapes, trial_shapes, test_elements, trial_elements, qd, quadstrat)
    verbose = (length(test_elements) > 256)
    myid = Threads.threadid()
    verbose && myid == 1 && print("dots out of 10: ")
    todo, done, pctg = length(test_elements), 0, 0
  for (p,(tcell,tptr)) in enumerate(zip(test_elements, test_cell_ptrs))
        for (q,(bcell,bptr)) in enumerate(zip(trial_elements, trial_cell_ptrs))

        fill!(zlocal, 0)
        # qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd, quadstrat)
        # momintegrals!(zlocal, biop,
        #     test_space,  tptr, tcell,
        #     trial_space, bptr, bcell, qrule)
        quadrule_momintegrals!(zlocal, biop,
            test_space,  tptr, tcell,
            trial_space, bptr, bcell, test_shapes, trial_shapes, p, q, qdl, quadstrat)
        I = length(test_assembly_data[p])
        J = length(trial_assembly_data[q])
        for j in 1 : J, i in 1 : I
            zij = zlocal[i,j]
            for (n,b) in trial_assembly_data[q][j]
                zb = zij*b
                for (m,a) in test_assembly_data[p][i]
                    store(a*zb, m, n)
        end end end end

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            verbose && myid == 1 && print(".")
            pctg = new_pctg
    end end
    verbose && myid == 1 && println("")
end
function assemblechunk_body_colored!(biop,
        test_space, testad, testelementids,
        trial_space, trialad, trialelementids,
        qd, store, scheduler, quadstrat::T) where {T<:RefQuadStrat}

    test_shapes = refspace(test_space)
    trial_shapes = refspace(trial_space)

    test_elements, test_assembly_data, test_cell_ptrs = testad
    trial_elements, trial_assembly_data, trial_cell_ptrs = trialad

    num_tshapes = numfunctions(test_shapes, domain(chart(geometry(test_space), first(geometry(test_space)))))
    num_bshapes = numfunctions(trial_shapes, domain(chart(geometry(trial_space), first(geometry(trial_space)))))

    @tasks for p in testelementids
        @set scheduler = scheduler
        @local begin zlocal = zeros(scalartype(biop, test_space, trial_space), num_tshapes, num_bshapes)
             qd = quaddatathreadlocal(biop, test_shapes, trial_shapes, test_elements, trial_elements, qd, quadstrat)
        end

        tcell, tptr = test_elements[p], test_cell_ptrs[p]
    
        for q in trialelementids
            bcell, bptr = trial_elements[q], trial_cell_ptrs[q]
            fill!(zlocal, 0)
            # @inline qrule = quadrule(biop, test_shapes, trial_shapes, p, tcell, q, bcell, qd, quadstrat)
            # momintegrals!(zlocal, biop,
            #     test_space,  tptr, tcell,
            #     trial_space, bptr, bcell, qrule)
            quadrule_momintegrals!(zlocal, biop,
                test_space,  tptr, tcell,
                trial_space, bptr, bcell, test_shapes, trial_shapes, p, q, qd, quadstrat)
            for j in 1:length(trial_assembly_data[q]), i in 1:length(test_assembly_data[p])
                zij = zlocal[i,j]
                for (n,b) in trial_assembly_data[q][j]
                    iszero(b) && continue
                    zb = zij*b
                    for (m,a) in test_assembly_data[p][i]
                        iszero(a) && continue
                        store(a*zb, m, n)
end end end end end end

# #### some dictionary checks 

# function moddict!(di, toadd)
#     di[:c] = toadd 
#     return nothing
# end
# function init_dict()
#     return Dict{Symbol,Int64}(:a => 1, :b => 2)
# end
# function run_dict(di, data)
#     out = 0
#     for dat in data 
#         moddict!(di, dat)
#         out += di[:c]
#     end
#     return out
# end

# #### some array checks 
# mutable struct SizeableArray{T}
#     data::Vector{T}
#     length::Int64
# end
# function Base.getindex(sa::SizeableArray, i)
#     return sa.data[i]
# end
# function Base.setindex!(sa::SizeableArray{T}, v, i) where {T}
#     if sa.length < i
#         push!(sa.data, zeros(T, i-sa.length)...)
#         sa.length = i
#     end
#     sa.data[i] = v
# end 

# ##### structured dictionary based on dispatch.
# struct StructDict{Param} end
# struct StructDictKey{Key} end
# mutable struct StructDictValue{Param, Key, Val} 
#     value::Val
# end
# Base.getindex