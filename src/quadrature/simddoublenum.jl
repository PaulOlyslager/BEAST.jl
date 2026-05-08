using LoopVectorization

abstract type Operations end
struct Cross <: Operations end
struct Dot <: Operations end
struct Times <: Operations end
struct STimesS <: Operations end
struct VTimesS <: Operations end
struct STimesV <: Operations end
function (::Cross)(a,b)
    return cross(a,b)
end
function (::Dot)(a,b)
    return transpose(a)*b
end
function (::Union{Times,STimesS,VTimesS,STimesV})(a,b) 
    return a*b
end
# struct SIMDDoubleNumQStrat{R} <: AbstractQuadStrat
#     outer_rule::R
#     inner_rule::R
#     # gre_cache::Matrix{Float64}
#     # gim_cache::Matrix{Float64}
#     # g_cache::Vector{Matrix{Float64}}
# end
# SIMDStratSupported = Union{DoubleNumQStrat, DoubleNumSauterQstrat, DoubleNumWiltonBogaertQStrat, DoubleNumWiltonSauterQStrat}
# struct SIMDQuadStrat{InnerStrat} <: AbstractQuadStrat
#     inner_strat::InnerStrat
#     function SIMDQuadStrat(inner_strat::T) where {T <: SIMDStratSupported}
#         return new(inner_strat)
#     end
#     function SIMDQuadStrat(inner_strat)
#         return inner_strat
#     end
# end



# function SIMDDoubleNumQStrat(outer_rule, inner_rule)
#     a = length(trgauss(outer_rule)[2])
#     b = length(trgauss(inner_rule)[2])
#     # gre_cache = zeros(a,b)
#     # gim_cache = zeros(a,b)
#     # g_cache = [zeros(a,b), zeros(a,b)]
#     SIMDDoubleNumQStrat(outer_rule, inner_rule)
# end
# function quaddata(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_elements, trial_elements, qs::SIMDDoubleNumQStrat)

#     # local_test_basis = refspace(test_basis)
#     # local_trial_basis = refspace(trial_basis)

#     test_quad_data  = quadpoints_simd(local_test_basis,  test_elements,  (qs.outer_rule,))
#     trial_quad_data = quadpoints_simd(local_trial_basis, trial_elements, (qs.inner_rule,))
#     # g_cache = generate_cacheSIMDDoubleNum(operator, test_quad_data, trial_quad_data)

#     return (tpoints = test_quad_data, bpoints = trial_quad_data)#, g_cache
# end
# function quaddata(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_elements, trial_elements, qs::SIMDQStrat)

#     # local_test_basis = refspace(test_basis)
#     # local_trial_basis = refspace(trial_basis)

#     qd = quaddata(operator, local_test_basis, local_trial_basis, test_elements, trial_elements, qs.inner_strat)
#     # g_cache = generate_cacheSIMDDoubleNum(operator, test_quad_data, trial_quad_data)

#     return (tpoints = qd.tpoints, bpoints = qd.bpoints, qd[2:end]...)#, g_cache
# end
# function quaddatathreadlocal(biop, tshapes, bshapes, test_elements, trial_elements, qd, quadstrat::SIMDDoubleNumQStrat)
#     g_cache = generate_cacheSIMDDoubleNum(biop, qd[1], qd[2])
#     return qd..., g_cache
# end
function generate_cacheSIMDDoubleNum(operator::CompDoubleKern{T,Op1,Op2,Kern}, test_quad_data, trial_quad_data) where {T,Op1,Kern,Op2}
    a = length(test_quad_data[1,1].weights)
    b = length(trial_quad_data[1,1].weights)
    return [zeros(a,b) for _ in 1:cachesizeSIMDDoubleNum(Kern)]
end
function generate_cacheSIMDDoubleNum(operator::MWSingleLayer3D, test_quad_data, trial_quad_data) 
    a = length(test_quad_data[1,1].weights)
    b = length(trial_quad_data[1,1].weights)
    return [zeros(a,b) for _ in 1:10]
end
cachesizeSIMDDoubleNum(::Type{<:HH3DGreen}) = 2
cachesizeSIMDDoubleNum(::Type{<:HH3DGradGreen}) = 6
function quadpoints_simd(f, els, rules, derivative_needed)
    if derivative_needed 
        return map(Iterators.product(rules, els)) do (rule, el)
        pws = quadpoints(el, rule)
        (weights = [w for (p,w) in pws], 
        points = [cartesian(p)[i] for i in 1:length(cartesian(pws[1][1])), (p,w) in pws], 
        values = [f(p)[j][1][i] for  i in 1:length(f(pws[1][1])[1][1]), (p,w) in pws,j in 1:length(f(pws[1][1]))], 
        derivatives =[f(p)[j][2][i] for  i in 1:length(f(pws[1][1])[1][2]), (p,w) in pws,j in 1:length(f(pws[1][1]))])
        end
    else
        return map(Iterators.product(rules, els)) do (rule, el)
            pws = quadpoints(el, rule)
            (weights = [w for (p,w) in pws], 
            points = [cartesian(p)[i] for i in 1:length(cartesian(pws[1][1])), (p,w) in pws], 
            values = [f(p)[j][1][i] for  i in 1:length(f(pws[1][1])[1][1]), (p,w) in pws,j in 1:length(f(pws[1][1]))])
        end
    end
    
end
# SIMDOperator(op::MWSingleLayer3D) = op
# SIMDOperator(op::CompDoubleKern) = op
# SIMDOperator(op) = error("SIMDOperator not implemented for operator of type $(typeof(op))")
requires_derivative(::MWSingleLayer3D) = true
requires_derivative(::CompDoubleKern) = false
# function quadrule(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,
#     quad_data, qs::SIMDQStrat{<:DoubleNumQStrat})

#     test_quad_rules  = quad_data[1]
#     trial_quad_rules = quad_data[2]

#     SIMDDoubleQuadRule(
#         test_quad_rules[1,test_id],
#         trial_quad_rules[1,trial_id],
#         SIMDOperator(operator),
#         # qs.gre_cache,
#         # qs.gim_cache,
#         quad_data[3]
#     )
# end

# function quadrule(operator::IntegralOperator,
#     local_test_basis, local_trial_basis,
#     test_id, test_element, trial_id, trial_element,
#     quad_data, qs::SIMDQStrat)
#     qe = quadrule(operator, local_test_basis, local_trial_basis, test_id, test_element, 
#         trial_id, trial_element, quad_data, qs.inner_strat)
#     SIMDRule(qe,quad_data,qs)
# end
# SIMDRulesSuported = Union{DoubleQuadRule}
# struct SIMDRule{InnerRule, BIOP, R}
#     innerrule::InnerRule
#     biop::BIOP 
#     g_cache::R
#     function SIMDRule(qr::T, qd, qs) where {T <: SIMDRulesSuported}
#         new(qr, SIMDOperator(biop), qd[end])
#     end
#     function SIMDRule(qr, qd, qs)
#         qr 
#     end
# end


struct SIMDDoubleQuadRule{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
#   biop::BIOP
  g_cache::Vector{Matrix{Float64}}
end
# SIMDDoubleQuadRule(outer_quad_points, inner_quad_points) = SIMDDoubleQuadRule(outer_quad_points, inner_quad_points, zeros(length(outer_quad_points.weights), length(inner_quad_points.weights)), zeros(length(outer_quad_points.weights), length(inner_quad_points.weights)))

"""
momintegrals!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
# function momintegrals!(biop::Maxwell3D.MWSingleLayer3D,
#     tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)

#     igd = Integrand(biop, tshs, bshs, tcell, bcell)

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points

#     # hoist split x/y/z arrays so indexing inside loops is only by wompid/wimpid
#     w_weights = womps.weights
#     w_pointsx = @view womps.points[1,:]
#     w_pointsy = @view womps.points[2,:]
#     w_pointsz = @view womps.points[3,:]
#     w_vals11 = @view womps.values[1,1,:]; w_vals12 = @view womps.values[1,2,:]; w_vals13 = @view womps.values[1,3,:]
#     w_vals21 = @view womps.values[2,1,:]; w_vals22 = @view womps.values[2,2,:]; w_vals23 = @view womps.values[2,3,:]
#     w_vals31 = @view womps.values[3,1,:]; w_vals32 = @view womps.values[3,2,:]; w_vals33 = @view womps.values[3,3,:]
#     w_der1 = @view womps.derivatives[1,:]; w_der2 = @view womps.derivatives[2,:]; w_der3 = @view womps.derivatives[3,:]

#     i_weights = wimps.weights
#     i_pointsx = @view wimps.points[1,:]
#     i_pointsy = @view wimps.points[2,:]
#     i_pointsz = @view wimps.points[3,:]
#     i_vals11 = @view wimps.values[1,1,:]; i_vals12 = @view wimps.values[1,2,:]; i_vals13 = @view wimps.values[1,3,:]
#     i_vals21 = @view wimps.values[2,1,:]; i_vals22 = @view wimps.values[2,2,:]; i_vals23 = @view wimps.values[2,3,:]
#     i_vals31 = @view wimps.values[3,1,:]; i_vals32 = @view wimps.values[3,2,:]; i_vals33 = @view wimps.values[3,3,:]
#     i_der1 = @view wimps.derivatives[1,:]; i_der2 = @view wimps.derivatives[2,:]; i_der3 = @view wimps.derivatives[3,:]

#     αim = imag(igd.operator.α)
#     αre = real(igd.operator.α)
#     βim = imag(igd.operator.β)
#     βre = real(igd.operator.β)
#     γim = imag(igd.operator.gamma)
#     γre = real(igd.operator.gamma)
#     z11re = 0.0
#     z12re = 0.0
#     z13re = 0.0

#     z21re = 0.0
#     z22re = 0.0
#     z23re = 0.0

#     z31re = 0.0
#     z32re = 0.0
#     z33re = 0.0

#     z11im = 0.0
#     z12im = 0.0
#     z13im = 0.0

#     z21im = 0.0
#     z22im = 0.0
#     z23im = 0.0

#     z31im = 0.0
#     z32im = 0.0
#     z33im = 0.0
#     @turbo for wompid in eachindex(w_weights)
#         tgeox = w_pointsx[wompid]
#         tgeoy = w_pointsy[wompid]
#         tgeoz = w_pointsz[wompid]
#         tvals1x = w_vals11[wompid]
#         tvals1y = w_vals12[wompid]
#         tvals1z = w_vals13[wompid]
#         tvals2x = w_vals21[wompid]
#         tvals2y = w_vals22[wompid]
#         tvals2z = w_vals23[wompid]
#         tvals3x = w_vals31[wompid]
#         tvals3y = w_vals32[wompid]
#         tvals3z = w_vals33[wompid]
#         tdivs1 = w_der1[wompid]
#         tdivs2 = w_der2[wompid]
#         tdivs3 = w_der3[wompid]
#         # M = length(tvals)
#         jx = w_weights[wompid]

#         for wimpid in eachindex(i_weights)
#             bgeox = i_pointsx[wimpid]
#             bgeoy = i_pointsy[wimpid]
#             bgeoz = i_pointsz[wimpid]
#             bvals1x = i_vals11[wimpid]
#             bvals1y = i_vals12[wimpid]
#             bvals1z = i_vals13[wimpid]
#             bvals2x = i_vals21[wimpid]
#             bvals2y = i_vals22[wimpid]
#             bvals2z = i_vals23[wimpid]
#             bvals3x = i_vals31[wimpid]
#             bvals3y = i_vals32[wimpid]
#             bvals3z = i_vals33[wimpid]
#             bdivs1 = i_der1[wimpid]
#             bdivs2 = i_der2[wimpid]
#             bdivs3 = i_der3[wimpid]
#             # N = length(bvals)
#             jy = i_weights[wimpid]

#             j = jx * jy
#             #### integrand evaluation


#             rx = tgeox-bgeox; ry = tgeoy-bgeoy; rz = tgeoz-bgeoz
#             R = sqrt(rx*rx + ry*ry + rz*rz)
#             iR = 1 / R
#             # green = exp(-γ*R)*(i4pi*iR)
#             greenre = cos(γim*R)*exp(-γre*R)*(i4pi*iR) 
#             greenim = - sin(γim*R)*exp(-γre*R)*(i4pi*iR)
#             αGre = αre * greenre - αim * greenim
#             αGim = αim * greenre + αre * greenim
#             βGre = βre * greenre - βim * greenim
#             βGim = βim * greenre + βre * greenim

#             # explicit 3x3 dot products between test vals and basis vals
#             d11 = tvals1x*bvals1x + tvals1y*bvals1y + tvals1z*bvals1z
#             d12 = tvals1x*bvals2x + tvals1y*bvals2y + tvals1z*bvals2z
#             d13 = tvals1x*bvals3x + tvals1y*bvals3y + tvals1z*bvals3z

#             d21 = tvals2x*bvals1x + tvals2y*bvals1y + tvals2z*bvals1z
#             d22 = tvals2x*bvals2x + tvals2y*bvals2y + tvals2z*bvals2z
#             d23 = tvals2x*bvals3x + tvals2y*bvals3y + tvals2z*bvals3z

#             d31 = tvals3x*bvals1x + tvals3y*bvals1y + tvals3z*bvals1z
#             d32 = tvals3x*bvals2x + tvals3y*bvals2y + tvals3z*bvals2z
#             d33 = tvals3x*bvals3x + tvals3y*bvals3y + tvals3z*bvals3z

#             z11re += j * (αGre * d11 + βGre * tdivs1 * bdivs1)
#             z12re += j * (αGre * d12 + βGre * tdivs1 * bdivs2)
#             z13re += j * (αGre * d13 + βGre * tdivs1 * bdivs3)

#             z11im += j * (αGim * d11 + βGim * tdivs1 * bdivs1)
#             z12im += j * (αGim * d12 + βGim * tdivs1 * bdivs2)
#             z13im += j * (αGim * d13 + βGim * tdivs1 * bdivs3)

#             z21re += j * (αGre * d21 + βGre * tdivs2 * bdivs1)
#             z22re += j * (αGre * d22 + βGre * tdivs2 * bdivs2)
#             z23re += j * (αGre * d23 + βGre * tdivs2 * bdivs3)

#             z21im += j * (αGim * d21 + βGim * tdivs2 * bdivs1)
#             z22im += j * (αGim * d22 + βGim * tdivs2 * bdivs2)
#             z23im += j * (αGim * d23 + βGim * tdivs2 * bdivs3)

#             z31re += j * (αGre * d31 + βGre * tdivs3 * bdivs1)
#             z32re += j * (αGre * d32 + βGre * tdivs3 * bdivs2)
#             z33re += j * (αGre * d33 + βGre * tdivs3 * bdivs3)

#             z31im += j * (αGim * d31 + βGim * tdivs3 * bdivs1)
#             z32im += j * (αGim * d32 + βGim * tdivs3 * bdivs2)
#             z33im += j * (αGim * d33 + βGim * tdivs3 * bdivs3)
#         end
#     end
#     z[1,1] += z11re + im*z11im
#     z[1,2] += z12re + im*z12im
#     z[1,3] += z13re + im*z13im
#     z[2,1] += z21re + im*z21im
#     z[2,2] += z22re + im*z22im
#     z[2,3] += z23re + im*z23im
#     z[3,1] += z31re + im*z31im
#     z[3,2] += z32re + im*z32im
#     z[3,3] += z33re + im*z33im
#     return z
# end

function momintegrals!(biop::MWSingleLayer3D,
    tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)

    igd = Integrand(biop, tshs, bshs, tcell, bcell)

    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points

    αim = imag(igd.operator.α)
    αre = real(igd.operator.α)
    βim = imag(igd.operator.β)
    βre = real(igd.operator.β)
    γim = imag(igd.operator.gamma)
    γre = real(igd.operator.gamma)
    # g_cache = (zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}), zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}))
    @turbo for wimpid in eachindex(wimps.weights)
            bgeoxi = wimps.points[1,wimpid]
            bgeoyi = wimps.points[2,wimpid]
            bgeozi = wimps.points[3,wimpid]
           for wompid in eachindex(womps.weights)
                tgeoxi = womps.points[1,wompid]
                tgeoyi = womps.points[2,wompid]
                tgeozi = womps.points[3,wompid]
                rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
                R = sqrt(rx*rx + ry*ry + rz*rz)
                iR = 1 / R
                e = exp(-γre*R)
                s, c = sincos(γim*R)
                scale = e * (i4pi * iR)
                greenre = c * scale
                greenim = -s * scale
                # green = exp(-γ*R)*(i4pi*iR)
                # strat.gre_cache[wompid, wimpid] = greenre
                # strat.gim_cache[wompid, wimpid] = greenim
                strat.g_cache[1][wompid, wimpid] = greenre 
                strat.g_cache[2][wompid, wimpid] = greenim
            end
        end
    for l in 1:3 
        for k in 1:3
            zre = 0.0
            zim = 0.0
            @turbo for wimpid in eachindex(wimps.weights)
                    bvalsxi = wimps.values[1,wimpid,l]
                    bvalsyi = wimps.values[2,wimpid,l]
                    bvalszi = wimps.values[3,wimpid,l]
                    bdiv = wimps.derivatives[wimpid,l]
                    jy = wimps.weights[wimpid]
            
            for wompid in eachindex(womps.weights)
                tvalsxi = womps.values[1,wompid,k]
                tvalsyi = womps.values[2,wompid,k]
                tvalszi = womps.values[3,wompid,k]
                tdiv = womps.derivatives[wompid,k]
                jx = womps.weights[wompid]
                    j = jx * jy
                    # greenre = strat.gre_cache[wompid, wimpid]
                    # greenim = strat.gim_cache[wompid, wimpid]
                    greenre = strat.g_cache[1][wompid, wimpid]
                    greenim = strat.g_cache[2][wompid, wimpid]
                    dotprod = tvalsxi*bvalsxi + tvalsyi*bvalsyi + tvalszi*bvalszi
                    tdivbdiv = tdiv * bdiv
                    A_re = αre*dotprod + βre*tdivbdiv
                    A_im = αim*dotprod + βim*tdivbdiv

                    zre += j * (greenre*A_re - greenim*A_im)
                    zim += j * (greenre*A_im + greenim*A_re)
                end
            end
 
            z[k,l] += zre + im*zim


        end
    end
   
    return z
end

TransposedStrat(a::SIMDDoubleQuadRule) = a

# function momintegrals!(biop::MWSingleLayer3D,
#     tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)

#     igd = Integrand(biop, tshs, bshs, tcell, bcell)

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points

#     αim = imag(igd.operator.α)
#     αre = real(igd.operator.α)
#     βim = imag(igd.operator.β)
#     βre = real(igd.operator.β)
#     γim = imag(igd.operator.gamma)
#     γre = real(igd.operator.gamma)
#     # g_cache = (zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}), zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}))
    
#     lwimps = length(wimps.weights)
#     lwomps = length(womps.weights)
#     ltot = lwimps * lwomps
#     for i in 0:(ltot ÷ 8)-1
#         bgeoxi = select_1djump8hor((@view wimps.points[1,:]), i, lwimps)
#         bgeoyi = select_1djump8hor((@view wimps.points[2,:]), i, lwimps)
#         bgeozi = select_1djump8hor((@view wimps.points[3,:]), i, lwimps)

#         tgeoxi = select_1djump8ver((@view womps.points[1,:]), i, lwimps)
#         tgeoyi = select_1djump8ver((@view womps.points[2,:]), i, lwimps)
#         tgeozi = select_1djump8ver((@view womps.points[3,:]), i, lwimps)
#         rx = tgeoxi - bgeoxi; ry = tgeoyi - bgeoyi; rz = tgeozi - bgeozi
#         R = sqrt(rx*rx + ry*ry + rz*rz)
#         iR = 1 / R
#         e = exp(-γre * R)
#         s, c = sin(γim * R), cos(γim * R)
#         scale = e * (i4pi * iR)
#         greenre = c * scale
#         greenim = -s * scale
#         # green = exp(-γ*R)*(i4pi*iR)
#         # strat.gre_cache[wompid, wimpid] = greenre
#         # strat.gim_cache[wompid, wimpid] = greenim
#         for l in 1:3
#             for k in 1:3
#                 bvalsxi = select_1djump8hor((@view wimps.values[1,:,l]), i, lwimps)
#                 bvalsyi = select_1djump8hor((@view wimps.values[2,:,l]), i, lwimps)
#                 bvalszi = select_1djump8hor((@view wimps.values[3,:,l]), i, lwimps)
#                 bdiv = select_1djump8hor((@view wimps.derivatives[1,:,l]), i, lwimps)
#                 jy = select_1djump8hor( wimps.weights, i, lwimps)

#                 tvalsxi = select_1djump8ver((@view womps.values[1,:,k]), i, lwimps)
#                 tvalsyi = select_1djump8ver((@view womps.values[2,:,k]), i, lwimps)
#                 tvalszi = select_1djump8ver((@view womps.values[3,:,k]), i, lwimps)
#                 tdiv = select_1djump8ver((@view womps.derivatives[1,:,k]), i, lwimps)
#                 jx = select_1djump8ver( womps.weights, i, lwimps)
#                 j = jx * jy
#                 # greenre = strat.gre_cache[wompid, wimpid]
#                 # greenim = strat.gim_cache[wompid, wimpid]
#                 # greenre = select_1djump8ver((@view strat.g_cache[1,:,:]), i, lwimps)
#                 # greenim = select_1djump8ver((@view strat.g_cache[2,:,:]), i, lwimps)
#                 dotprod = tvalsxi*bvalsxi + tvalsyi*bvalsyi + tvalszi*bvalszi
#                 tdivbdiv = tdiv * bdiv
#                 A_re = αre*dotprod + βre*tdivbdiv
#                 A_im = αim*dotprod + βim*tdivbdiv
#                 zre = sum(j * (greenre*A_re - greenim*A_im))
#                 zim = sum(j * (greenre*A_im + greenim*A_re))
#                 z[k,l] += zre + im*zim
#             end
#         end

#     end
#     return z 
# end
# function momintegrals!(biop::MWSingleLayer3D,
#     tshs, bshs, tcell, bcell, z, strat::SIMDDoubleQuadRule)

#     igd = Integrand(biop, tshs, bshs, tcell, bcell)

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points

#     αim = imag(igd.operator.α)
#     αre = real(igd.operator.α)
#     βim = imag(igd.operator.β)
#     βre = real(igd.operator.β)
#     γim = imag(igd.operator.gamma)
#     γre = real(igd.operator.gamma)
#     # g_cache = (zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}), zero(MMatrix{length(womps.weights), length(wimps.weights), Float64}))
    
#     lwimps = length(wimps.weights)
#     lwomps = length(womps.weights)
#     ltot = lwimps * lwomps
#     # assign everything in cache first in matrix form.
#     for i in eachindex(womps.weights)
#         for j in eachindex(wimps.weights)
#             strat.g_cache[1][i,j] = wimps.points[1,i]
#             strat.g_cache[2][i,j] = wimps.points[2,i]
#             strat.g_cache[3][i,j] = wimps.points[3,i]
#             strat.g_cache[4][i,j] = womps.points[1,j]
#             strat.g_cache[5][i,j] = womps.points[2,j]
#             strat.g_cache[6][i,j] = womps.points[3,j]
#             strat.g

#     for i in 0:(ltot ÷ 8)-1
#         bgeoxi = select_1djump8hor((@view wimps.points[1,:]), i, lwimps)
#         bgeoyi = select_1djump8hor((@view wimps.points[2,:]), i, lwimps)
#         bgeozi = select_1djump8hor((@view wimps.points[3,:]), i, lwimps)

#         tgeoxi = select_1djump8ver((@view womps.points[1,:]), i, lwimps)
#         tgeoyi = select_1djump8ver((@view womps.points[2,:]), i, lwimps)
#         tgeozi = select_1djump8ver((@view womps.points[3,:]), i, lwimps)
#         rx = tgeoxi - bgeoxi; ry = tgeoyi - bgeoyi; rz = tgeozi - bgeozi
#         R = sqrt(rx*rx + ry*ry + rz*rz)
#         iR = 1 / R
#         e = exp(-γre * R)
#         s, c = sin(γim * R), cos(γim * R)
#         scale = e * (i4pi * iR)
#         greenre = c * scale
#         greenim = -s * scale
#         # green = exp(-γ*R)*(i4pi*iR)
#         # strat.gre_cache[wompid, wimpid] = greenre
#         # strat.gim_cache[wompid, wimpid] = greenim
#         for l in 1:3
#             for k in 1:3
#                 bvalsxi = select_1djump8hor((@view wimps.values[1,:,l]), i, lwimps)
#                 bvalsyi = select_1djump8hor((@view wimps.values[2,:,l]), i, lwimps)
#                 bvalszi = select_1djump8hor((@view wimps.values[3,:,l]), i, lwimps)
#                 bdiv = select_1djump8hor((@view wimps.derivatives[1,:,l]), i, lwimps)
#                 jy = select_1djump8hor( wimps.weights, i, lwimps)

#                 tvalsxi = select_1djump8ver((@view womps.values[1,:,k]), i, lwimps)
#                 tvalsyi = select_1djump8ver((@view womps.values[2,:,k]), i, lwimps)
#                 tvalszi = select_1djump8ver((@view womps.values[3,:,k]), i, lwimps)
#                 tdiv = select_1djump8ver((@view womps.derivatives[1,:,k]), i, lwimps)
#                 jx = select_1djump8ver( womps.weights, i, lwimps)
#                 j = jx * jy
#                 # greenre = strat.gre_cache[wompid, wimpid]
#                 # greenim = strat.gim_cache[wompid, wimpid]
#                 # greenre = select_1djump8ver((@view strat.g_cache[1,:,:]), i, lwimps)
#                 # greenim = select_1djump8ver((@view strat.g_cache[2,:,:]), i, lwimps)
#                 dotprod = tvalsxi*bvalsxi + tvalsyi*bvalsyi + tvalszi*bvalszi
#                 tdivbdiv = tdiv * bdiv
#                 A_re = αre*dotprod + βre*tdivbdiv
#                 A_im = αim*dotprod + βim*tdivbdiv
#                 zre = sum(j * (greenre*A_re - greenim*A_im))
#                 zim = sum(j * (greenre*A_im + greenim*A_re))
#                 z[k,l] += zre + im*zim
#             end
#         end

#     end
#     return z 
# end
#     # for wimpid in eachindex(wimps.weights)
#     #         bgeoxi = wimps.points[1,wimpid]
#     #         bgeoyi = wimps.points[2,wimpid]
#     #         bgeozi = wimps.points[3,wimpid]
#     #        for wompid in eachindex(womps.weights)
#     #             tgeoxi = womps.points[1,wompid]
#     #             tgeoyi = womps.points[2,wompid]
#     #             tgeozi = womps.points[3,wompid]
#     #             rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
#     #             R = sqrt(rx*rx + ry*ry + rz*rz)
#     #             iR = 1 / R
#     #             e = exp(-γre*R)
#     #             s, c = sincos(γim*R)
#     #             scale = e * (i4pi * iR)
#     #             greenre = c * scale
#     #             greenim = -s * scale
#     #             # green = exp(-γ*R)*(i4pi*iR)
#     #             # strat.gre_cache[wompid, wimpid] = greenre
#     #             # strat.gim_cache[wompid, wimpid] = greenim
#     #             strat.g_cache[1][wompid, wimpid] = greenre 
#     #             strat.g_cache[2][wompid, wimpid] = greenim
#     #         end
#     #     end
#     # for l in 1:3 
#     #     for k in 1:3
#     #         zre = 0.0
#     #         zim = 0.0
#     #         @turbo for wimpid in eachindex(wimps.weights)
#     #                 bvalsxi = wimps.values[1,wimpid,l]
#     #                 bvalsyi = wimps.values[2,wimpid,l]
#     #                 bvalszi = wimps.values[3,wimpid,l]
#     #                 bdiv = wimps.derivatives[wimpid,l]
#     #                 jy = wimps.weights[wimpid]
            
#     #         for wompid in eachindex(womps.weights)
#     #             tvalsxi = womps.values[1,wompid,k]
#     #             tvalsyi = womps.values[2,wompid,k]
#     #             tvalszi = womps.values[3,wompid,k]
#     #             tdiv = womps.derivatives[wompid,k]
#     #             jx = womps.weights[wompid]
#     #                 j = jx * jy
#     #                 # greenre = strat.gre_cache[wompid, wimpid]
#     #                 # greenim = strat.gim_cache[wompid, wimpid]
#     #                 greenre = strat.g_cache[1][wompid, wimpid]
#     #                 greenim = strat.g_cache[2][wompid, wimpid]
#     #                 dotprod = tvalsxi*bvalsxi + tvalsyi*bvalsyi + tvalszi*bvalszi
#     #                 tdivbdiv = tdiv * bdiv
#     #                 A_re = αre*dotprod + βre*tdivbdiv
#     #                 A_im = αim*dotprod + βim*tdivbdiv

#     #                 zre += j * (greenre*A_re - greenim*A_im)
#     #                 zim += j * (greenre*A_im + greenim*A_re)
#     #             end
#     #         end
 
#     #         z[k,l] += zre + im*zim


#     #     end
#     # end
   
# #     return z
# # end
# using SIMD

@generated function select_1djump8ver(a::AbstractVector{T},  index, ly) where {T}
    ex = quote
    end

    # build expressions for the 8 elements of each tuple
    ai_elems = Expr[:(VecElement(a[((8*index + $k) ÷ ly) + 1])) for k in 1:8]
    # bi_elems = Expr[:(b[mod1(8*index + $k, ly)]) for k in 1:8]

    ai_tuple = Expr(:tuple, ai_elems...)
    # bi_tuple = Expr(:tuple, bi_elems...)

    push!(ex.args, :(ai = $ai_tuple))
    # push!(ex.args, :(bi = $bi_tuple))
    push!(ex.args, quote return Vec(ai) end)

    return ex
end
@generated function select_1djump8hor( b::AbstractVector{T}, index, ly) where {T}
    ex = quote
    end

    # build expressions for the 8 elements of each tuple
    # ai_elems = Expr[:(a[((8*index + $k) ÷ ly) + 1]) for k in 1:8]
    bi_elems = Expr[:(VecElement(b[mod1(8*index + $k, ly)])) for k in 1:8]

    # ai_tuple = Expr(:tuple, ai_elems...)
    bi_tuple = Expr(:tuple, bi_elems...)

    # push!(ex.args, :(ai = $ai_tuple))
    push!(ex.args, :(bi = $bi_tuple))
    push!(ex.args, quote return Vec(bi) end)

    return ex
end
#### SIMD for composed operators:

#### create a green-function cache based on the green function.
#### evaluate the green function in the points and fill the cache (inline so use generated function)
#### loop over test and trial elements then loop over points and fill z_local, apply operations. for type operations I can implement generated functions 
#### but then user needs to use those type generated functions.

@generated function momintegrals!(biop::CompDoubleKern{T,Op1,Op2,Kern},
    tshs::TSpace, bshs::BSpace, tcell, bcell, z, strat::SIMDDoubleQuadRule) where {T,Op1,Kern,Op2,TSpace,BSpace}
    M = numfunctions(tshs, domain(tcell))
    N = numfunctions(bshs, domain(bcell))
    Core.println("M = $M, N = $N")
    # load variables
    vars = _load_variables(Kern)
    # green function evaluation loop with cache filling
    kernexp = _kernexp(Kern)
    # expressions of the two operators 
    ## possible operators: cross, dot, times (scalar*scalar), times (scalar*vector), times (vector*scalar)
    ## green function is put in lx, ly, lz  or l 
    ## basis function is put in rx, ry, rz, of r 
    ## maps to solx, soly, solz, or sol 
    assign_green_realexp = assign_green_real(Kern)
    assign_green_imagexp = assign_green_imag(Kern)
    get_greenrealexp = get_greenreal(Kern)
    get_greenimagexp = get_greenimag(Kern)
    op1 = operationexp(Op1,Kern,bshs)
    op2 = operationexp(Op2,tshs,Kern)
    trialspacetovarmapexp = trialspacetovarmap(BSpace)
    testspacetovarmapexp = testspacetovarmap(TSpace)
    loadvaluestestexp = loadvalues_test(TSpace)
    loadvaluestrialexp = loadvalues_trial(BSpace)
    ex = quote 
        womps = strat.outer_quad_points
        wimps = strat.inner_quad_points
        $vars
       @turbo for wimpid in eachindex(wimps.weights)
            bgeoxi = wimps.points[1,wimpid]
            bgeoyi = wimps.points[2,wimpid]
            bgeozi = wimps.points[3,wimpid]
           for wompid in eachindex(womps.weights)
                tgeoxi = womps.points[1,wompid]
                tgeoyi = womps.points[2,wompid]
                tgeozi = womps.points[3,wompid]

                $kernexp
                $assign_green_realexp
                $assign_green_imagexp

           end
        end
        for li in 1:$N 
            for ki in 1:$M
                zre = 0.0
                zim = 0.0
            @turbo for wimpid in eachindex(wimps.weights)
                            $loadvaluestrialexp
                            # bdiv = wimps.derivatives[wimpid,li]
                            jy = wimps.weights[wimpid]
                    
                            for wompid in eachindex(womps.weights)
                                $loadvaluestestexp
                                # tdiv = womps.derivatives[wompid,ki]
                                jx = womps.weights[wompid]
                                j = jx * jy

                                $get_greenrealexp
                                $trialspacetovarmapexp
                                $op2
                                $testspacetovarmapexp
                                $op1
                                zre += j * r

                                $get_greenimagexp
                                $trialspacetovarmapexp
                                $op2
                                $testspacetovarmapexp
                                $op1
                                zim += j * r

                            end
                        end
                z[ki,li] += zre + im*zim
            end

        end

            return z
    end
    # display(ex)
    # error("stop")
    return ex
    # zlocal filling

    # this is the same as the above function but with the operator application split into two steps and with the operator parameters hoisted out of the loops. I can also implement generated functions for the operator application to get more type stability but then users need to use those generated functions. 
end
valuetype(::Type{RTRefSpace{T}}) where {T} = SVector{3,T}
CompScienceMeshes.domain(a::Type{CompScienceMeshes.Simplex{A,B,C,D,T}}) where {A,B,C,D,T} = CompScienceMeshes.ReferenceSimplex{B} 
function _load_variables(::Type{U}) where {U <: Union{HH3DGreen, HH3DGradGreen}}
    quote
        γim = imag(biop.kernel.gamma)
        γre = real(biop.kernel.gamma)
    end
end
function loadvalues_test(::Type{U}) where {U}
    if valuetype(U) <: Number
        return quote
            tvals = womps.values[1, wompid,ki]
        end
    end
    if valuetype(U) <: SVector{3}
        return quote
            tvalsxi = womps.values[1,wompid,ki]
            tvalsyi = womps.values[2,wompid,ki]
            tvalszi = womps.values[3,wompid,ki]
        end
    end
end
function loadvalues_trial(::Type{U}) where {U}
    if valuetype(U) <: Number
        return quote
            bvals = wimps.values[1,wimpid,li]
        end
    end
    if valuetype(U) <: SVector{3}
        return quote
            bvalsxi = wimps.values[1,wimpid,li]
            bvalsyi = wimps.values[2,wimpid,li]
            bvalszi = wimps.values[3,wimpid,li]
        end
    end
end
function _kernexp(::Type{<:HH3DGreen})
    quote
        rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
        R = sqrt(rx*rx + ry*ry + rz*rz)
        iR = 1 / R
        e = exp(-γre*R)
        s, c = sincos(γim*R)
        scale = e * (i4pi * iR)
        greenre = c * scale
        greenim = -s * scale
    end
end
function _kernexp(::Type{<:HH3DGradGreen})
    quote
        rx = tgeoxi-bgeoxi; ry = tgeoyi-bgeoyi; rz = tgeozi-bgeozi
        R = sqrt(rx*rx + ry*ry + rz*rz)
        iR = 1 / R
        e = exp(-γre*R)
        s, c = sincos(γim*R)
        scale = e * (i4pi * iR)
        greenre = c*scale
        greenim = -s*scale

        gradgreenre = iR*((-γre-iR)*greenre + γim*greenim)
        gradgreenim = iR*((-γre-iR)*greenim - γim*greenre)

       
        greenrex = gradgreenre * rx
        greenrey = gradgreenre * ry
        greenrez = gradgreenre * rz
        greenimx = gradgreenim * rx
        greenimy = gradgreenim * ry
        greenimz = gradgreenim * rz
    end
end
operationexp(::Type{Cross}, Kern, bshs) = _cross_exp()
operationexp(::Type{Dot}, Kern, bshs) = _dot_exp()
operationexp(::Type{STimesS}, Kern, bshs) = _scalar_times_scalar_exp()
operationexp(::Type{VTimesS}, Kern, bshs) = _vector_times_scalar_exp()
operationexp(::Type{STimesV}, Kern, bshs) = _scalar_times_vector_exp()

function trialspacetovarmap(a::Type{T}) where {T <: RefSpace}

    if _size(valuetype(T)) == (3,)
        ex = quote
            rx = bvalsxi; ry = bvalsyi; rz = bvalszi
        end
    elseif _size(valuetype(T)) == 1
        ex = quote
            r = bvals
        end
    else
        # k = _size(valuetype(T))
        # Core.println(k)
        error("unsupported trial space value type")
    end
    return ex
end
function testspacetovarmap(a::Type{T}) where {T <: RefSpace}

    if _size(valuetype(T)) == (3,)
        ex = quote
            lx = tvalsxi; ly = tvalsyi; lz = tvalszi
        end
    elseif _size(valuetype(T)) == 1
        ex = quote
            l = tvals
        end
    else
        error("unsupported trial space value type")
    end
    return ex
end
valuetype(::Type{_LocalBasisDot{T}}) where {T} = T
function valuetype(::Type{_LocalBasisTimes{T,U,V}}) where {T,U,V}

    valuetype(U) <: Number && valuetype(V) <: Number && return promote_type(U,V)
    return SVector{3,T}
end
valuetype(::Type{_LocalBasisCross{T,U,V}}) where {T,U,V} = SVector{3,T}
valuetype(::Type{NormalVector}) = nothing
# function operationexp(::Type{Times}, Kern, bshs)
#     left_vector = _size(returntype(Kern)) == 3
#     right_vector = _size(returntype(bshs)) == 3
# end
_size(a::Type{<:SVector}) = size(a)
_size(a::Type{<:Number}) = 1
function _scalar_times_scalar_exp()
    quote
        sol = l*r
        r = sol
    end
end
function _scalar_times_vector_exp()
    quote
        solx = l*rx
        soly = l*ry
        solz = l*rz
        rx = solx
        ry = soly
        rz = solz
    end
end
function _vector_times_scalar_exp()
    quote
        solx = lx*r
        soly = ly*r
        solz = lz*r
        rx = solx
        ry = soly
        rz = solz
    end
end
function _dot_exp()
    quote
        r = rx*lx + ry*ly + rz*lz
    end
end
function _cross_exp()
    quote
        solx = ry*lz - rz*ly
        soly = rz*lx - rx*lz
        solz = rx*ly - ry*lx
        rx = solx
        ry = soly
        rz = solz
    end
end

function get_greenreal(kern::Type{<:HH3DGreen})
    quote
        l = strat.g_cache[1][wompid, wimpid]
    end
end
function get_greenimag(kern::Type{<:HH3DGreen})
    quote
        l = strat.g_cache[2][wompid, wimpid]
    end
end
function get_greenreal(kern::Type{<:HH3DGradGreen})
    quote
        lx = strat.g_cache[1][wompid, wimpid]
        ly = strat.g_cache[2][wompid, wimpid]
        lz = strat.g_cache[3][wompid, wimpid]
    end
end
function get_greenimag(kern::Type{<:HH3DGradGreen})
    quote
        lx = strat.g_cache[4][wompid, wimpid]
        ly = strat.g_cache[5][wompid, wimpid]
        lz = strat.g_cache[6][wompid, wimpid]
    end 
end
function assign_green_real(kern::Type{<:HH3DGreen})
    quote
        strat.g_cache[1][wompid, wimpid] = greenre
    end
end
function assign_green_imag(kern::Type{<:HH3DGreen})
    quote
        strat.g_cache[2][wompid, wimpid] = greenim
    end
end
function assign_green_real(kern::Type{<:HH3DGradGreen})
    quote
        strat.g_cache[1][wompid, wimpid] = greenrex
        strat.g_cache[2][wompid, wimpid] = greenrey
        strat.g_cache[3][wompid, wimpid] = greenrez
    end
end
function assign_green_imag(kern::Type{<:HH3DGradGreen})
    quote
        strat.g_cache[4][wompid, wimpid] = greenimx
        strat.g_cache[5][wompid, wimpid] = greenimy
        strat.g_cache[6][wompid, wimpid] = greenimz
    end
end