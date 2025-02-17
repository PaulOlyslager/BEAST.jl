# T: coeff type
# Degree: degree
# Dim1: dimension of the support + 1
mutable struct LagrangeRefSpace{T,Degree,Dim1,NF} <: RefSpace{T,NF} end

numfunctions(s::LagrangeRefSpace{T,D,2}) where {T,D} = D+1
numfunctions(s::LagrangeRefSpace{T,0,3}) where {T} = 1
numfunctions(s::LagrangeRefSpace{T,1,3}) where {T} = 3
numfunctions(s::LagrangeRefSpace{T,2,3}) where {T} = 6

valuetype(ref::LagrangeRefSpace{T}, charttype) where {T} =
        SVector{numfunctions(ref), Tuple{T,T}}

# Evaluate constant lagrange elements on anything
(ϕ::LagrangeRefSpace{T,0})(tp) where {T} = SVector(((value=one(T), derivative=zero(T)),))

# Evaluate linear Lagrange elements on a segment
function (f::LagrangeRefSpace{T,1,2})(mp) where T
    u = mp.bary[1]
    j = jacobian(mp)
    SVector(
        (value=  u, derivative=-1/j),
        (value=1-u, derivative= 1/j))
end

# Evaluete linear lagrange elements on a triangle
function (f::LagrangeRefSpace{T,1,3})(t) where T
    u,v,w, = barycentric(t)
    # SVector(
    #     (value=u,),
    #     (value=v,),
    #     (value=w,))

    j = jacobian(t)
    p = t.patch
    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
    SVector(
        (value=u, curl=σ*(p[3]-p[2])/j),
        (value=v, curl=σ*(p[1]-p[3])/j),
        (value=w, curl=σ*(p[2]-p[1])/j))
end


"""
    f(tangent_space, Val{:withcurl})

Compute the values of the shape functions together with their curl.
"""
# function (f::LagrangeRefSpace{T,1,3})(t, ::Type{Val{:withcurl}}) where T
#     # Evaluete linear Lagrange elements on a triange, together with their curl
#     j = jacobian(t)
#     u,v,w, = barycentric(t)
#     p = t.patch
#     σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
#     SVector(
#         (value=u, curl=σ*(p[3]-p[2])/j),
#         (value=v, curl=σ*(p[1]-p[3])/j),
#         (value=w, curl=σ*(p[2]-p[1])/j)
#     )
# end


# Evaluate constant Lagrange elements on a triangle, with their curls
function (f::LagrangeRefSpace{T,0,3})(t, ::Type{Val{:withcurl}}) where T
    i = one(T)
    z = zero(cartesian(t))
    SVector(((value=i, curl=z,),))
end


function curl(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), -sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), +sh.coeff)
    return [sh1, sh2]
end

function gradient(ref::LagrangeRefSpace{T,1,4}, sh, tet) where {T}
    this_vert = tet.vertices[sh.refid]
    # other_verts = deleteat(tet.vertices, sh.refid)
    # opp_face = simplex(other_verts...)
    opp_face = faces(tet)[sh.refid]
    ctr_opp_face = center(opp_face)
    n = normal(ctr_opp_face)
   # n = normalize(tangents(ctr_opp_face,1)×tangents(ctr_opp_face,2))
    h = -dot(this_vert - cartesian(ctr_opp_face), n)
    @assert h > 0
    gradval = -(1/h)*n
    output = Vector{Shape{T}}()
    for (i,edge) in enumerate(CompScienceMeshes.edges(tet))
        ctr_edge = center(edge)
        tgt = tangents(ctr_edge,1)
        tgt = normalize(tgt)
        lgt = volume(edge)
        cff = -lgt * dot(tgt, gradval)
        isapprox(cff, 0, atol=sqrt(eps(T))) && continue
        push!(output, Shape(sh.cellid, i, sh.coeff * cff))
    end
    return output
end

function gradient(ref::LagrangeRefSpace{T,1,3} where {T}, sh, el)
    sh1 = Shape(sh.cellid, mod1(sh.refid+1,3), +sh.coeff)
    sh2 = Shape(sh.cellid, mod1(sh.refid+2,3), -sh.coeff)
    return [sh1, sh2]
end

function gradient(ref::LagrangeRefSpace{T,1,2}, sh, seg) where {T}

    sh.refid == 1 && return [Shape(sh.cellid, 1, +sh.coeff/volume(seg))]
    @assert sh.refid == 2
    return [Shape(sh.cellid, 1, -sh.coeff/volume(seg))]

end



function strace(x::LagrangeRefSpace, cell, localid, face)

    Q = zeros(scalartype(x),2,3)

    p1 = neighborhood(face, 1)
    p2 = neighborhood(face, 0)

    u1 = carttobary(cell, cartesian(p1))
    u2 = carttobary(cell, cartesian(p2))

    P1 = neighborhood(cell, u1)
    P2 = neighborhood(cell, u2)

    vals1 = x(P1)
    vals2 = x(P2)

    for j in 1:numfunctions(x)
        Q[1,j] = vals1[j].value
        Q[2,j] = vals2[j].value
    end

    Q
end


function restrict(refs::LagrangeRefSpace{T,0}, dom1, dom2) where T
    #Q = eye(T, numfunctions(refs))
    Q = Matrix{T}(I, numfunctions(refs), numfunctions(refs))
end

function restrict(f::LagrangeRefSpace{T,1}, dom1, dom2) where T

    D = numfunctions(f)
    Q = zeros(T, D, D)

    # for each point of the new domain
    for i in 1:D
        v = verticeslist(dom2)[i]

        # find the barycentric coordinates in dom1
        uvn = carttobary(dom1, v)

        # evaluate the shape functions in this point
        x = neighborhood(dom1, uvn)
        fx = f(x)

        for j in 1:D
            Q[j,i] = fx[j][1]
        end
    end

    return Q
end



## Quadratic Lagrange element on a triangle
function (f::LagrangeRefSpace{T,2,3})(t) where T
    j = jacobian(t)
    p = t.patch

    #curl=(p[3]-p[2])/j),
     #   (value=v, curl=(p[1]-p[3])/j),
      #  (value=w, curl=(p[2]-p[1])/j)

    σ = sign(dot(normal(t), cross(p[1]-p[3],p[2]-p[3])))
     SVector(
        (value=u*(2*u-1), curl=σ*(p[3]-p[2])*(4u-1)/j),
        (value=v*(2*v-1), curl=σ*(p[1]-p[3])*(4v-1)/j),
        (value=w*(2*w-1), curl=σ*(p[2]-p[1])*(4w-1)/j),
        (value=4*v*w, curl=4*σ*(w*(p[1]-p[3])+v*(p[2]-p[1]))/j),
        (value=4*w*u, curl=4*σ*(w*(p[3]-p[2])+u*(p[2]-p[1]))/j),
        (value=4*u*v, curl=4*σ*(u*(p[1]-p[3])+v*(p[3]-p[2]))/j),
    )
end

# function (f::LagrangeRefSpace{T,2,3})(t, ::Type{Val{:withcurl}}) where T
#     # Evaluete quadratic Lagrange elements on a triange, together with their curl
#     tu = tangents(t,1)
#     tv = tangents(t,2)
#     j = jacobian(t)*sign(dot(normal(t),tu×tv))
#     u,v,w, = barycentric(t)
#     p = t.patch
#     SVector(
#         (value=u*(2*u-1), curl=(p[3]-p[2])*(4u-1)/j),
#         (value=v*(2*v-1), curl=(p[1]-p[3])*(4v-1)/j),
#         (value=w*(2*w-1), curl=(p[2]-p[1])*(4w-1)/j),
#         (value=4*v*w, curl=4*(w*(p[1]-p[3])+v*(p[2]-p[1]))/j),
#         (value=4*w*u, curl=4*(w*(p[3]-p[2])+u*(p[2]-p[1]))/j),
#         (value=4*u*v, curl=4*(u*(p[1]-p[3])+v*(p[3]-p[2]))/j),
#     )
# end

function curl(ref::LagrangeRefSpace{T,2,3} where {T}, sh, el)
    #curl of lagc0d2 as combination of bdm functions 
    z=zero(typeof(sh.coeff))
    if sh.refid < 4
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+1,6), +sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+2,6), -3*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+3,6), +3*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+4,6), -sh.coeff)
    else
        sh1 = Shape(sh.cellid, mod1(2*sh.refid+4,6), z*sh.coeff)
        sh2 = Shape(sh.cellid, mod1(2*sh.refid+5,6), -4*sh.coeff)
        sh3 = Shape(sh.cellid, mod1(2*sh.refid+6,6), +4*sh.coeff)
        sh4 = Shape(sh.cellid, mod1(2*sh.refid+7,6), z*sh.coeff)
    end
    return [sh1, sh2, sh3, sh4]
end


const _vert_perms_lag = [
    (1,2,3),
    (2,3,1),
    (3,1,2),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

const _dof_perms_lag0 = [
    (1),
    (1),
    (1),
    (1),
    (1),
    (1),
]
const _dof_perms_lag1 = [
    (1,2,3),
    (3,1,2),
    (2,3,1),
    (2,1,3),
    (1,3,2),
    (3,2,1),
]

function dof_permutation(::LagrangeRefSpace{<:Any,0}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag0[i]
end

function dof_permutation(::LagrangeRefSpace{<:Any,1}, vert_permutation)
    i = findfirst(==(tuple(vert_permutation...)), _vert_perms_lag)
    return _dof_perms_lag1[i]
end