struct NonConformingTouchQRule{S}
    conforming_qstrat::S
    test_overlapping_edge_index::Int
    bsis_overlapping_edge_index::Int
end

function momintegrals!(op, test_locspace, bsis_locspace,
        τ::CompScienceMeshes.Simplex, σ::CompScienceMeshes.Simplex,
        out, qrule::NonConformingTouchQRule)

    T = coordtype(τ)
    P = eltype(τ.vertices)

    i = qrule.test_overlapping_edge_index
    j = qrule.bsis_overlapping_edge_index

    τs, σs = _conforming_refinement_touching_triangles(τ,σ,i,j)

    # test conformity
    for a in τs
        for b in σs
            @assert _test_conformity(a, b)
        end
    end

    @assert volume(τ) ≈ sum(volume.(τs))
    @assert volume(σ) ≈ sum(volume.(σs))
 
    qstrat = qrule.conforming_qstrat
    qdata = quaddata(op, test_locspace, bsis_locspace, τs, σs, qstrat)

    for (p,tchart) in enumerate(τs)
        for (q,bchart) in enumerate(σs)
            qrule = quadrule(op, test_locspace, bsis_locspace,
                p, tchart, q, bchart, qdata, qstrat)
            # @show qrule

            P = restrict(test_locspace, τ, tchart)
            Q = restrict(bsis_locspace, σ, bchart)
            zlocal = zero(out)

            momintegrals!(op, test_locspace, bsis_locspace,
                tchart, bchart, zlocal, qrule)

            # out .+= P * zlocal * Q'
            for i in axes(P,1)
                for j in axes(Q,1)
                    for k in axes(P,2)
                        for l in axes(Q,2)
                            out[i,j] += P[i,k] * zlocal[k,l] * Q[j,l]
end end end end end end end


function _conforming_refinement_touching_triangles(τ,σ,i,j)
    λ = faces(τ)[i]
    μ = faces(σ)[j]
    ρ = CompScienceMeshes.intersection(λ,μ)[1]

    τ_verts = [τ[mod1(i+2,3)], τ[mod1(i,3)], τ[mod1(i+1,3)]]
    σ_verts = [σ[mod1(j+2,3)], σ[mod1(j,3)], σ[mod1(j+1,3)]]

    T = coordtype(τ)
    P = eltype(τ.vertices)

    U = T[]
    V = P[]
    for v in ρ.vertices
        if CompScienceMeshes.isinside(λ,v)
            push!(V,v)
            push!(U, carttobary(λ,v)[1])
    end end
    if length(U) == 2
        if U[1] < U[2]
            temp = V[1]
            V[1] = V[2]
            V[2] = temp
    end end
    append!(τ_verts, V)

    U = T[]
    V = P[]
    for v in ρ.vertices
        if CompScienceMeshes.isinside(μ,v)
            push!(V,v)
            push!(U, carttobary(μ,v)[1])
    end end
    if length(U) == 2
        if U[1] < U[2]
            temp = V[1]
            V[1] = V[2]
            V[2] = temp
    end end
    append!(σ_verts, V)

    τ_verts = push!(τ_verts[2:end], τ_verts[1])
    σ_verts = push!(σ_verts[2:end], σ_verts[1])

    # @show τ_verts

    τ_charts = [ simplex(τ_verts[1], τ_verts[i], τ_verts[i+1]) for i in 2:length(τ_verts)-1 ]
    σ_charts = [ simplex(σ_verts[1], σ_verts[i], σ_verts[i+1]) for i in 2:length(σ_verts)-1 ]

    signs = Int.(sign.(dot.(normal.(τ_charts),Ref(normal(τ)))))
    τ_charts = flip_normal.(τ_charts,signs)
    signs = Int.(sign.(dot.(normal.(σ_charts),Ref(normal(σ)))))
    σ_charts = flip_normal.(σ_charts,signs)

    return τ_charts, σ_charts
end


function _num_common_vertices(τ, σ)
    hits = 0
    T = coordtype(σ)
    tol = eps(T) * 10^3
    for v in τ.vertices
        for w in σ.vertices
            if norm(v - w) < tol
                hits += 1
                break
    end end end
    return hits
end


function _test_conformity(τ, σ)
    if CompScienceMeshes.overlap(τ, σ)
        # if _num_common_vertices(τ, σ) != 3
        #     @infiltrate
        # end
        return _num_common_vertices(τ, σ) == 3
    end

    for eτ in faces(τ)
        for eσ in faces(σ)
            if CompScienceMeshes.overlap(eτ, eσ)
                # if _num_common_vertices(τ, σ) != 2
                #     @infiltrate
                # end
                return _num_common_vertices(τ, σ) == 2
            end
    end end

    for u in τ.vertices
        for v in σ.vertices
            su = simplex([u], Val{0})
            sv = simplex([v], Val{0})
            if CompScienceMeshes.overlap(su, sv)
                # if _num_common_vertices(τ, σ) != 1
                #     @infiltrate
                # end
                return _num_common_vertices(τ, σ) == 1
    end end end

    @assert _num_common_vertices(τ, σ) == 0
    return _num_common_vertices(τ, σ) == 0
end