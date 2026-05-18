##### implementation of the rules of:
# M. M. Botha, "A Family of Augmented Duffy Transformations for Near-Singularity 
# Cancellation Quadrature," in IEEE Transactions on Antennas and Propagation, 
# vol. 61, no. 6, pp. 3123-3134, June 2013, doi: 10.1109/TAP.2013.2252137.


using FastGaussQuadrature
import CompScienceMeshes.jacobian
function project_point(s::CompScienceMeshes.Simplex, point)
    verticess = CompScienceMeshes.vertices(s)
    t1 = normalize(verticess[2] - verticess[1])
    t2 = normalize(verticess[3] - verticess[1])
    t2p = normalize(t2 - dot(t2, t1) * t1)
    a1 = dot(point-verticess[1], t1)
    a2 = dot(point-verticess[1], t2p)
    vp = a1 * t1 + a2 * t2p + verticess[1]
    return vp
end



struct BothaTriangle
    ex  
    ey
    v1
    xrange 
    yrange
    sign
end
function surface(bt::BothaTriangle)
    return (bt.xrange[2]-bt.xrange[1]) * (bt.yrange[2]-bt.yrange[1]) / 2 * bt.sign
end
function BothaTriangle(vertices, normal)
    ex = normalize(vertices[2] - vertices[3])
    ey = normalize((vertices[3] - vertices[1]) + dot(vertices[1] - vertices[3], ex) * ex)
    v1 = vertices[1]
    xrange = (dot(vertices[3] - vertices[1], ex), dot(vertices[2] - vertices[1], ex))
    @assert xrange[1] < xrange[2]
    yrange = (0.0, dot(vertices[3] - vertices[1], ey))
    @assert yrange[1] < yrange[2] yrange
    return BothaTriangle(ex, ey, v1, xrange, yrange, sign(dot(cross(vertices[2] - vertices[1], vertices[3] - vertices[1]), normal)))
end
function cartesian_coords(bt::BothaTriangle, u, v)
    return bt.v1 + u * bt.ex + v * bt.ey
end

abstract type BothaRule end 

struct ADR1C <: BothaRule end
function coordmap(::ADR1C)
    return [x-> asinh(x), (y,w,z) -> sqrt(y^2 + (z/cosh(w))^2)]
end
function invmap(::ADR1C)
    return  [(w,v,z) -> sqrt(v^2-(z/cosh(w))^2)*sinh(w), (w,v,z) -> sqrt(v^2-(z/cosh(w))^2)]
end
function jacobian(::ADR1C)
    return (x,y,z,v,w) -> sqrt(x^2+y^2+z^2)
end
struct ADR1L <: BothaRule end
function coordmap(::ADR1L)
    return [x -> asinh(x), (y,w,z) -> asinh(y*cosh(w)/z)]
end
function invmap(::ADR1L)
    return [(w,v,z) -> sinh(w)*sinh(v)*z/cosh(w), (w,v,z) -> sinh(v)*z/cosh(w)]
end
function jacobian(::ADR1L)
    return (x,y,z,v,w) -> sqrt(x^2 + y^2 + z^2) * z *sinh(v)/cosh(w)
end
struct ADR2C <: BothaRule end
function coordmap(::ADR2C)
    return [x -> atan(x), (y,w,z) -> log(sqrt(y^2 + (z*cos(w))^2))]
end
function invmap(::ADR2C)
    return [(w,v,z) -> sqrt(exp(2*v)-(z*cos(w))^2)*tan(w), (w,v,z) -> sqrt(exp(2*v)-(z*cos(w))^2)]
end
function jacobian(::ADR2C)
    return (x,y,z,v,w) -> (x^2 + y^2 + z^2)
end
struct ADR2L <: BothaRule end
function coordmap(::ADR2L)
    return [x -> asinh(x), (y,w,z) -> cosh(w)/z*atan(y*cosh(w)/z)]
end
function invmap(::ADR2L)
    return [(w,v,z) -> sinh(w)*tan(v*z/cosh(w))*z/cosh(w), (w,v,z) -> tan(v*z/cosh(w))*z/cosh(w)]
end
function jacobian(::ADR2L)
    return (x,y,z,v,w) -> (x^2 + y^2 + z^2) * z * tan(v*z/cosh(w))/cosh(w)^2 
end
struct ADR3C <: BothaRule end
function coordmap(::ADR3C)
    return [x -> atan(x), (y,w,z) -> -1/sqrt(y^2 + (z*cos(w))^2)]
end
function invmap(::ADR3C)
    return [(w,v,z) -> sqrt(1/v^2-(z*cos(w))^2)*tan(w), (w,v,z) -> sqrt(1/v^2-(z*cos(w))^2)]
end
function jacobian(::ADR3C)
    return (x,y,z,v,w) -> (x^2 + y^2 + z^2)^(3/2) * cos(w)
end
struct ADR3L <: BothaRule end
function coordmap(::ADR3L)
    return [x -> asinh(x), (y,w,z) -> (cosh(w)/z)^2 * y / sqrt(y^2 + (z/cosh(w))^2)]
end
function invmap(::ADR3L)
    return [(w,v,z) -> (z/cosh(w))/sqrt((cosh(w)/z)^4/v^2-1)*sinh(w), (w,v,z) -> (z/cosh(w))/sqrt((cosh(w)/z)^4/v^2-1)]
end
function jacobian(::ADR3L)
    return (x,y,z,v,w) -> (x^2 + y^2 + z^2)^(3/2) * z * v / cosh(w)^3 / sqrt((cosh(w)/z)^4 - v^2)
end
struct BothaTriangleQuadStrat <: AbstractQuadStrat
    botharule:: BothaRule
    innerrule1
    innerrule2
end
quaddata(op,rs,els,qs::BothaTriangleQuadStrat) = nothing
function quadrule(op,refspace,p,y,q,el,qdata,qs::BothaTriangleQuadStrat)
    ### create Botha triangles
    bts = subdivide_triangle(el, y)
    ### set up quadrature
    p1,w1 = FastGaussQuadrature.gausslegendre(qs.innerrule1)
    p2,w2 = FastGaussQuadrature.gausslegendre(qs.innerrule2)
    ptot, wtot, vals = typeof(center(el))[], typeof(w1[1])[], typeof(refspace(center(el)))[]
    for bt in bts
        x1,x2 = bt.xrange
        y1,y2 = bt.yrange
        z = norm(cartesian(y)-bt.v1 - dot(cartesian(y)-bt.v1, bt.ex) * bt.ex - dot(cartesian(y)-bt.v1, bt.ey) * bt.ey)
        p1new, b1new = rescale_to_interval(p1, w1, coordmap(qs.botharule)[1](x1/y2), coordmap(qs.botharule)[1](x2/y2))
        for (w,weight1) in zip(p1new, b1new)
            p2new, b2new = rescale_to_interval(p2, w2, coordmap(qs.botharule)[2](0, w,z), coordmap(qs.botharule)[2](y2, w, z))
            for (v, weight2) in zip(p2new, b2new)
                xi,yi = invmap(qs.botharule)[1](w,v,z), invmap(qs.botharule)[2](w,v,z)
                point = cartesian_coords(bt, xi, yi)
                weight = jacobian(qs.botharule)(xi,yi,z,v,w) * bt.sign * weight1 * weight2
                push!(ptot, neighborhood(el,carttobary(el,point)))
                push!(wtot, weight)
                push!(vals, refspace(ptot[end]))

            end
        end
    end
  #  @assert sum(wtot) ≈ sum(surface.(bts)) (sum(wtot) , sum(surface.(bts)))
    return [(value = v, point = p, weight = w) for (v,p,w) in zip(vals, ptot, wtot)]
end

function rescale_to_interval(points, weights, a, b)
    newpoints = ((b - a) / 2) * points .+ (a + b) / 2
    newweights = ((b - a) / 2) * weights
    return newpoints, newweights
end

function subdivide_triangle(simplex, point)
    vp = project_point(simplex, point)
    n = normal(simplex)
    v1, v2, v3 = vertices(simplex)
    bts = [BothaTriangle([vi, vj, vk], n) for (vi, vj, vk) in
           [(vp, v1, v2),
           (vp, v2, v3),
           (vp, v3, v1)] if norm(cross(vj - vi, vk - vi)) > 1e-14]
    return bts
end