using BEAST
using CompScienceMeshes
using StaticArrays
using Plots
using LinearAlgebra

Γ = Mesh([(@SVector [0.0, 0.0, 0.0]),
          (@SVector [1.0, 0.0, 0.0]),
          (@SVector [0.0, 1.0, 0.0])],
         [(@SVector [1, 2, 3])])
X = BEAST.lagrangecxd0(Γ)

G0 = BEAST.HH3DGreen(0.0)
∇G0 = BEAST.HH3DGradGreen(0.0)
u = BEAST.PotentialIntegralOperator{2}(G0,*,B->B)
xrange = [0.5]
yrange = range(-1.0, stop=1.0, length=20)
zrange = [0.001]
pts = [point(x,y,z) for x ∈ xrange, y ∈ yrange, z ∈ zrange];

n = 8

qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR1C(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR1C", title="Potential along line x=0.5, z=$(zrange[1])")
qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR1L(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot!(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR1L", title="Potential along line x=0.5, z=$(zrange[1])")
qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR2C(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot!(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR2C", title="Potential along line x=0.5, z=$(zrange[1])")
qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR2L(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot!(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR2L", title="Potential along line x=0.5, z=$(zrange[1])")
qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR3C(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot!(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR3C", title="Potential along line x=0.5, z=$(zrange[1])")
qs = BEAST.BothaTriangleQuadStrat(
BEAST.ADR3L(),
    n,
    n)
vals = potential(u, pts, [1.0], X; type=Float64, quadstrat = qs)
plot!(yrange, vals[1,:,1]; xlabel="y", ylabel="Potential", label="ADR3L", title="Potential along line x=0.5, z=$(zrange[1])")





vals = potential(u, pts, [1.0], X;  quadstrat = qs)





using LegendrePolynomials
using Plots
using FastGaussQuadrature
function plot_legendre(l)
    trange = range(-1.0, stop=1.0, length=10000)
    vals = Pl.(trange,l)
    display(plot(trange, vals))
end

function fit_polynomial(f, deg, points)
    A = zeros(length(points), deg+1)
    for (i,p) in enumerate(points)
        for j in 0:deg
            A[i,j+1] = p^j
        end
    end
    b = [f(p) for p in points]
    coeffs = A\b
    return x -> sum(c*x^(j-1) for (j,c) in enumerate(coeffs))
end

function cont_function(datapoints, xpoints)
    function f(x)
        for i in 1:length(xpoints)-1
            if xpoints[i] ≤ x ≤ xpoints[i+1]
                t = (x - xpoints[i]) / (xpoints[i+1] - xpoints[i])
                return (1-t)*datapoints[i] + t*datapoints[i+1]
            end
        end
        error("x out of bounds")
    end
    return f
end

#### without cut
l = 30
cvals = cont_function(vals, yrange)
p,w = FastGaussQuadrature.gausslegendre(l)
fitted = fit_polynomial(cvals, l-1, p)
plot(trange,fitted.(trange))
plot!(trange,cvals.(trange))

#### with cut
l = 10
c1 = 0.0
c2 = 0.5
p1,w1 = BEAST.rescale_to_interval(p,w,-1.0,c1)
p2,w2 = BEAST.rescale_to_interval(p,w,c1, c2)
p3,w3 = BEAST.rescale_to_interval(p,w,c2, 1.0)
fitted1 = fit_polynomial(cvals, l-1, p1)
fitted2 = fit_polynomial(cvals, l-1, p2)
fitted3 = fit_polynomial(cvals, l-1, p3)
trange1 = range(-1.0, stop=c1, length=1000)
trange2 = range(c1, stop=c2, length=1000)
trange3 = range(c2, stop=1.0, length=1000)
plot(trange1,fitted1.(trange1))
plot!(trange2,fitted2.(trange2))
plot!(trange3,fitted3.(trange3))
plot!(trange,cvals.(trange))