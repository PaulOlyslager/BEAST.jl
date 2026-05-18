using CompScienceMeshes, BEAST, LinearAlgebra
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = meshsphere(radius=1.0, h=0.25)

X = raviartthomas(Γ)

Δt = 0.1
Nt = 30
T = timebasisc0d1(Δt, Nt)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration = 2 * 20 * Δt
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)

@hilbertspace j
@hilbertspace j′

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)

tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
@time A = assemble(SL, W, V)

T = eltype(A)
S = zeros(T, size(A)[1:2])
BEAST.ConvolutionOperators.timeslice!(S, A, 1)

iS = inv(S)
b = assemble(E, W)

nt = numfunctions(temporalbasis(V))
@time marchonintime(iS, A, b, nt)


@profview xefie = BEAST.motsolve(tdefie)

import PlotlyJS
PlotlyJS.plot(xefie[1,:])
fcr, geo = facecurrents(xefie[:,125], X)
PlotlyJS.plot(patch(geo, norm.(fcr)))


function test()
    A = rand(1000,1000)
    v =rand(1000)
    Binv = lu(A)
    for i in 1:200
        v = Binv\v
        for j in 1:21
            v = A*v
        end
    end
end

@time test()

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)

fcr, geo = facecurrents(ue, X)
PlotlyJS.plot(patch(geo, norm.(fcr)))