struct DoubleQuadRule{P,Q}
  outer_quad_points::P
  inner_quad_points::Q
end


"""
momintegrals!(biop, tshs, bshs, tcell, bcell, interactions, strat)

Function for the computation of moment integrals using simple double quadrature.
"""
function momintegrals!(biop,
    tshs, bshs, tcell, bcell, z, strat::DoubleQuadRule)

    igd = Integrand(biop, tshs, bshs, tcell, bcell)

    womps = strat.outer_quad_points
    wimps = strat.inner_quad_points
    
    for womp in womps
        tgeo = womp.point
        tvals = womp.value
        M = length(tvals)
        jx = womp.weight
        
        for wimp in wimps
            bgeo = wimp.point
            bvals = wimp.value
            N = length(bvals)
            jy = wimp.weight

            j = jx * jy

            z1 = j * igd(tgeo, bgeo, tvals, bvals)
            for n in 1:N
                for m in 1:M
                    z[m,n] += z1[m,n]
            end end
        end
    end

    return z
 end
# struct Batch{T,N}
#     data::SVector{N,T}
# end

# function momintegrals!(biop,
#     tshs, bshs, tcell, bcell, z, strat::DoubleQuadRuleSimd)

#     igd = Integrand(biop, tshs, bshs, tcell, bcell)

#     womps = strat.outer_quad_points
#     wimps = strat.inner_quad_points
#     numwimps = length(wimps)
#     numbatches = numwimps ÷ 8 + 1
#     #create_batches
#     wimp_batches = [Batch(SVector{8}(wimps[i:min(i+7,numwimps)]...)) for i in 1:8:numwimps]

#     for womp in womps
#         tgeo = womp.point
#         tvals = womp.value
#         M = length(tvals)
#         jx = womp.weight
        
#         for wimpbatch in wimp_batches
#             bgeo_batch = wimp.point
#             bvals_batch = wimp.value
#             N = length(bvals)
#             jy = wimp.weight

#             j = jx * jy

#             z1 = j * igd(tgeo, bgeo, tvals, bvals)
#             for n in 1:N
#                 for m in 1:M
#                     z[m,n] += z1[m,n]
#             end end
#         end
#     end

#     return z
# end

# _TransposedStrat(a::DoubleQuadRule) = a


# function f1(a::Vector{SVector{3,ComplexF64}}, b::Vector{SVector{3,ComplexF64}})
#     out = zero(ComplexF64)
#     for i in 1:length(a)
#         out += dot(a[i], b[i])
#     end
#     return out
# end
# function f2(a::Vector{SVector{3,ComplexF64}}, b::Vector{SVector{3,ComplexF64}})
#     out = zero(ComplexF64)
#     numbatches = length(a) ÷ 8
#     for i in 0:numbatches-1
#         out += sum(dot.((@view a[i*8+1:i*8+8]), (@view b[i*8+1:i*8+8])))
#     end
#     return out
# end