mutable struct NDLCCBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end
redefine_geometrie(basis::NDLCCBasis{T,M,P},geo::F) where {T,M,P,F}= NDLCCBasis{T,F,P}(geo,basis.fns,basis.pos)

NDLCCBasis(geo, fns) = NDLCCBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(space::NDLCCBasis{T}) where {T} = NDLCCRefSpace{T}()
Base.similar(space::NDLCCBasis{T,M,P} where {T,M,P}, geo, fns, pos) = NDLCCBasis(geo, fns, pos)

"""
nedelecc3d(mesh, edges::Array{SArray{Tuple{2},Int64,1,2},1})

returns an object of the type NDLCCBasis.
"""

function nedelecc3d(mesh, edges)
    T = coordtype(mesh)
    P = vertextype(mesh)
    num_edges = numcells(edges)

    C = connectivity(edges, mesh, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    fns = Vector{Vector{Shape{T}}}(undef,num_edges)
    pos = Vector{P}(undef,num_edges)
    for (i,e) in enumerate(edges)

        fns[i] = Vector{Shape{T}}()
        pos[i] = cartesian(center(chart(edges,e)))

        for k in nzrange(C,i)

            j = rows[k]
            s = vals[k]
            push!(fns[i], Shape{T}(j, abs(s), sign(s)))
        end
    end

    NDLCCBasis(mesh, fns, pos)
end

"""
nedelecc3d(mesh) returns an object of type NDLCCBasis
"""

function nedelecc3d(mesh)
    edges = skeleton(mesh,1)
    nedelecc3d(mesh, edges)
end


function nedelec(mesh::CompScienceMeshes.AbstractMesh{U,4} where {U},
    edges=skeleton(mesh,1))

    nedelecc3d(mesh, edges)
end

curl(space::NDLCCBasis, geo, fns) = NDLCDBasis(geo, fns, space.pos)

ttrace(X::NDLCCBasis, geo, fns) = RTBasis(geo, fns, deepcopy(X.pos))
