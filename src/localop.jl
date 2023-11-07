####### designing neighborhood object which can save 2 normals
# using ReusePatterns
# using CompScienceMeshes
# struct MeshPointNormals
#     meshpoint::CompScienceMeshes.MeshPointNM
#     testnormal
#     trialnormal
# end
#ReusePatterns.@forward (MeshPointNormals, :meshpoint) CompScienceMeshes.MeshPointNM

# code = forward((MeshPointNormals, :meshpoint),CompScienceMeshes.MeshPointNM)
# uit = [code[1]]
# for line in code[2:length(code)]
#     push!(uit,replace(line, " CompScienceMeshes " => " "))
# end
# eval.(Meta.parse.(uit))

# function extendedneighborhood(p::CompScienceMeshes.Simplex, bary,nt,nb)
#     MeshPointNormals(CompScienceMeshes.neighborhood(p,bary),nt,nb)
# end

# normals(m::MeshPointNormals) = (m.testnormal,m.trialnormal)
using CollisionDetection


abstract type LocalOperator <: Operator end


function allocatestorage(op::LocalOperator, test_functions, trial_functions,
    storage_trait::Type{Val{:bandedstorage}}, longdelays_trait)

    T = scalartype(op, test_functions, trial_functions)

    M = Int[]
    N = Int[]
    V = T[]

    function store(v,m,n)
        push!(M,m)
        push!(N,n)
        push!(V,v)
    end

    function freeze()
        nrows = numfunctions(test_functions)
        ncols = numfunctions(trial_functions)
        return sparse(M,N,V, nrows, ncols)
    end

    return freeze, store
end

function allocatestorage(op::LocalOperator, testfunctions, trialfunctions,
    storage_trait::Type{Val{:sparsedicts}}, longdelays_trait)

    T = scalartype(op, testfunctions, trialfunctions)

    m = numfunctions(testfunctions)
    n = numfunctions(trialfunctions)
    Z = SparseMatrixDict{T,Int}(m,n)

    store(v,m,n) = (Z[m,n] += T(v))
    freeze() = SparseArrays.SparseMatrixCSC(Z)

    return freeze, store
end

function allocatestorage(op::LocalOperator, test_functions, trial_functions,
    storage_trait::Type{Val{:densestorage}}, longdelays_trait)

    T = scalartype(op, test_functions, trial_functions)

    Z = zeros(T, numfunctions(test_functions), numfunctions(trial_functions))
    store(v,m,n) = (Z[m,n] += T(v))
    freeze() = Z

    return freeze, store
end

function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store,
        threading::Type{Threading{:multi}};
        quadstrat=defaultquadstrat(biop, tfs, bfs))

    if geometry(tfs) == geometry(bfs)
        @warn "Function cannot be used for infinitely thin surfaces!!!"
        return assemble_local_matched!(biop, tfs, bfs, store; quadstrat)
    end

    if CompScienceMeshes.refines(geometry(tfs), geometry(bfs))
        return assemble_local_refines!(biop, tfs, bfs, store; quadstrat)
    end

    return assemble_local_mixed!(biop, tfs, bfs, store; quadstrat)
end

function assemble!(biop::LocalOperator, tfs::Space, bfs::Space, store,
    threading::Type{Threading{:single}};
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    if geometry(tfs) == geometry(bfs)
        @warn "Function cannot be used for infinitely thin surfaces!!!"
        return assemble_local_matched!(biop, tfs, bfs, store; quadstrat)
    end

    if CompScienceMeshes.refines(geometry(tfs), geometry(bfs))
        return assemble_local_refines!(biop, tfs, bfs, store; quadstrat)
    end

    return assemble_local_mixed!(biop, tfs, bfs, store; quadstrat)
end

function assemble_local_matched!(biop::LocalOperator, tfs::Space, bfs::Space, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    verbose = length(tels) > 10_000
    verbose && print(string(typeof(biop))*" dots out of 20: ")
    todo, done, pctg = length(tels), 0, 0
   # locmat = zeros(scalartype(biop, trefs, brefs), numfunctions(trefs), numfunctions(brefs))
    for (p,tcell) in enumerate(tels)

        P = ta2g[p]
        bcell = bels[P]
        
        @assert same_cell(bcell,tcell)
        q = bg2a[P]
        q == 0 && continue

        qr = quadrule(biop, trefs, brefs, tcell, bcell,qd, quadstrat)
  #      fill!(locmat, 0)
        locmat = cellinteractions(biop, trefs, brefs, tcell, bcell, qr)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                
                store(a * locmat[i,j] * b, m, n)
        
        end end

        new_pctg = round(Int, (done += 1) / todo * 100)
        verbose && new_pctg > pctg + 4 && (print("."); pctg = new_pctg)
    end
end


function assemble_local_refines!(biop::LocalOperator, tfs::Space, bfs::Space, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    println("Using 'refines' algorithm for local assembly:")

    tgeo = geometry(tfs)
    bgeo = geometry(bfs)
    @assert CompScienceMeshes.refines(tgeo, bgeo)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tels, tad, ta2g = assemblydata(tfs)
    bels, bad, ba2g = assemblydata(bfs)

    bg2a = zeros(Int, length(geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    print(string(typeof(biop))*" dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        P = ta2g[p]
        Q = CompScienceMeshes.parent(tgeo, P)
        q = bg2a[Q]

        bcell = bels[q]
        @assert overlap(tcell, bcell)
        #s = sign(dot(normal(tcell),normal(bcell)))
        isct = intersection2(tcell, bcell)
        for (tcell2,bcell2) in isct

            P = restrict(brefs, bcell, bcell2)
            Q = restrict(trefs, tcell, tcell2)

            qr = quadrule(biop, trefs, brefs, tcell2,bcell2, qd, quadstrat)
            zlocal = cellinteractions(biop, trefs, brefs, tcell2, bcell2, qr)
            zlocal = Q * zlocal * P'

            for i in 1 : numfunctions(trefs)
                for j in 1 : numfunctions(brefs)
                    for (m,a) in tad[p,i]
                        for (n,b) in bad[q,j]
                            store(a * zlocal[i,j] * b, m, n)
                        end # next basis function this cell supports
                    end # next test function this cell supports
                end # next refshape on basis side
            end # next refshape on test side

        end # next cell in intersection

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")

end

function assemble_local_matched!(biop::LocalOperator, tfs::subdBasis, bfs::subdBasis, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs))

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)
    for (p,tcell) in enumerate(tels)
        bcell = bels[p]
        @assert same_cell(bcell,tcell)
        qr = quadrule(biop, trefs, brefs, tcell,bcell, qd, quadstrat)
        locmat = cellinteractions(biop, trefs, brefs, tcell,bcell, qr)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p][i], (n,b) in bad[p][j]
                store(a * locmat[i,j] * b, m, n)
end end end end
function same_cell(a,b)
sort(verticeslist(a)) ≈ sort(verticeslist(b))
end

function elementstree(elements)

    nverts = dimension(eltype(elements)) + 1
    ncells = length(elements)

    @assert !isempty(elements)
    P = eltype(verticeslist(elements[1]))
    T = coordtype(eltype(elements))

    points = zeros(P, ncells)
    radii = zeros(T, ncells)

    for i in 1 : ncells

        verts = verticeslist(elements[i])

        bary = verts[1]
        for j in 2:length(verts)
            bary += verts[j]
        end

        points[i] = bary / nverts
        for j in 1 : nverts
            radii[i] = max(radii[i], norm(verts[j]-points[i]))
        end
    end

    return Octree(points, radii)
end


"""
    assemble_local_mixed(biop::LocalOperator, tfs, bfs)

For use when basis and test functions are defined on different meshes
"""
function assemble_local_mixed!(biop::LocalOperator, tfs::Space{T}, bfs::Space{T}, store;
    quadstrat=defaultquadstrat(biop, tfs, bfs)) where {T}

    tol = sqrt(eps(T))

    trefs = refspace(tfs)
    brefs = refspace(bfs)

    tels, tad = assemblydata(tfs)
    bels, bad = assemblydata(bfs)

    qd = quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    # store the bcells in an octree
    tree = elementstree(bels)

    print(string(typeof(biop))*" dots out of 10: ")
    todo, done, pctg = length(tels), 0, 0
    for (p,tcell) in enumerate(tels)

        tc, ts = boundingbox(verticeslist(tcell))
        pred = (c,s) -> boxesoverlap(c,s,tc,ts)

        for box in boxes(tree, pred)
            for q in box
                bcell = bels[q]

                if overlap(tcell, bcell)

                    isct = intersection2(tcell, bcell)
                    for (cellt,cellb) in isct
                        volume(cellt) < tol && continue

                        P = restrict(brefs, bcell, cellb)
                        Q = restrict(trefs, tcell, cellt)

                        qr = quadrule(biop, trefs, brefs, cellt,cellb, qd, quadstrat)
                        zlocal = cellinteractions(biop, trefs, brefs, cellt, cellb, qr)
                        zlocal = Q * zlocal * P'

                        for i in 1 : numfunctions(trefs)
                            for j in 1 : numfunctions(brefs)
                                for (m,a) in tad[p,i]
                                    for (n,b) in bad[q,j]
                                        store(a * zlocal[i,j] * b, m, n)
                                    end # next basis function this cell supports
                                end # next test function this cell supports
                            end # next refshape on basis side
                        end # next refshape on test side

                    end # next cell in intersection
                end # if overlap
            end # next cell in the basis geometry
        end # next box in the octree

        done += 1
        new_pctg = round(Int, done / todo * 100)
        if new_pctg > pctg + 9
            print(".")
            pctg = new_pctg
        end
    end # next cell in the test geometry

    println("")
end


# function cellinteractions_matched!(zlocal, biop, trefs, brefs, cell, qr,tcell=nothing,bcell=nothing)

#     num_tshs = length(qr[1][3])
#     num_bshs = length(qr[1][4])

#     # zlocal = zeros(Float64, num_tshs, num_bshs)
#     for q in qr

#         w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
#         j = w * jacobian(mp)
#         kernel = kernelvals(biop, mp)
        
#         for n in 1 : num_bshs
#             bval = bvals[n]
#             for m in 1 : num_tshs
#                 tval = tvals[m]

#                 igd = integrand(biop, kernel, mp, tval, bval)
#                 zlocal[m,n] += j * igd
#             end
#         end
#     end

#     return zlocal
# end

function cellinteractions_matched!(zlocal, biop, trefs, brefs, cellt,cellb, qr)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for q in qr

        w, (tp,bp), tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * jacobian(tp)
        kernel = kernelvals(biop, tp)
        
        zlocal += j * integrand(biop, kernel, tp, bp, tvals, bvals)
    end
    return zlocal
end

# function cellinteractions(biop, trefs::U, brefs::V, cell,qr,tcell,bcell) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}
    
#     num_tshs = length(qr[1][3])
#     num_bshs = length(qr[1][4])

#     zlocal = zeros(T, num_tshs, num_bshs)
     
#     for q in qr

#         w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
#         j = w * jacobian(mp)
#         kernel = kernelvals(biop, mp)

#         for m in 1 : num_tshs
#             tval = tvals[m]

#             for n in 1 : num_bshs
#                 bval = bvals[n]

#                 igd = integrand(biop, kernel, mp, tval, bval)
#                 zlocal[m,n] += j * igd

#             end
#         end
#     end

#     return zlocal
# end

function cellinteractions(biop, trefs::U, brefs::V, cellt,cellb,qr) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}
    
    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    
    for q in qr


        w, (tp,bp), tvals, bvals = q[1], q[2], q[3], q[4]
        
        j = w * jacobian(tp)
        kernel = kernelvals(biop, tp)

        zlocal += j * integrand(biop, kernel, tp, bp, tvals, bvals)
        
    end

    return zlocal
end
function getvalue(list::Matrix{T}) where {T}
    
    return SVector{length(list),T}([i for i in list])
end
# function cellinteractions(biop, trefs::U, brefs::V, cell,qr) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}
    
#     num_tshs = length(qr[1][3])
#     num_bshs = length(qr[1][4])

#     zlocal = zeros(T, num_tshs, num_bshs)
     
#     for q in qr

#         w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
#         j = w * jacobian(mp)
#         kernel = kernelvals(biop, mp)

#         for m in 1 : num_tshs
#             tval = tvals[m]

#             for n in 1 : num_bshs
#                 bval = bvals[n]

#                 igd = integrand(biop, kernel, mp, tval, bval)
#                 zlocal[m,n] += j * igd

#             end
#         end
#     end

#     return zlocal
# end

