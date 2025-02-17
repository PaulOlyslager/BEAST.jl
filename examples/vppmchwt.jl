using BEAST
using CompScienceMeshes
using LinearAlgebra
using AbstractTrees
using StaticArrays
using Plots
function T1(κ;trace=true,tr=1,kwargs...)
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)

    end
    if trace
        int = -[γₛ(∇Gx,tr) -γₛ(Gn,tr) ;
        ZeroOperator() -τ(∇Gdotn,tr) ]
        
    else
        int = -[nt×(∇Gx) -nt×(Gn) ;
        ZeroOperator() -(∇Gdotn) ]
    end
    return BEAST.matrix_to_bilform(int)
end
function T2(κ;trace=true,tr=1,kwargs...)
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)
    end

    if trace
        int = -[γₛ(∇Gx,tr) ZeroOperator();
        γₙ(Gr,tr) -γₙ(∇G,tr)]
        
    else
        int = -[nt×(∇Gx) ZeroOperator();
        nt ⋅(Gr) -nt ⋅(∇G)]
    end
    return BEAST.matrix_to_bilform(int)
end
function A(κ;trace=true,tr=1,kwargs...)
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)

    end
    if trace
        int = -[γₛ(Gr,tr) -γₛ(∇G,tr);
        τ(∇Gdot,tr)  κ^2*τ(Gr,tr)]
        
    else
        int = -[nt×(Gr) -nt×(∇G);
        (∇Gdot)  κ^2*(Gr)]
    end
    return BEAST.matrix_to_bilform(int)
end

function Bb(κ;trace=true,tr=1,kwargs...)
    G = BEAST.greenhh3d(wavenumber=κ)
    gradG = BEAST.∇(G)
    if trace
        ∇Gx =  BEAST.build_potential(gradG×B)
        Gn = BEAST.build_potential(G*(n*B))
        #Gnx = BEAST.build_potential(G*(n × B))
        ∇G = BEAST.build_potential(gradG*B)
        ∇Gdotn = BEAST.build_potential(nb ⋅ (gradG*B))
        ∇Gdot = BEAST.build_potential(B ⋅ gradG)
        
        Gr = BEAST.build_potential(G*B)
        ∇G∇B = BEAST.build_potential(gradG*BEAST.div(B))
        ∇Gxn = BEAST.build_potential(gradG×(n*B))
    else
        ∇Gx =  gradG×B
        Gn = G*(nb*B)
        #Gnx = G*(n × B)
        ∇G = gradG*B
        ∇Gdotn = nb ⋅ (gradG*B)
        ∇Gdot = B ⋅ gradG
        
        Gr = G*B
        ∇G∇B = gradG*BEAST.div(B)
        ∇Gxn = gradG×(nb*B)

    end
    if trace
        int = -[γₛ(∇G∇B,tr)+κ^2*γₛ(Gr,tr) -γₛ(∇Gxn,tr);
        γₙ(∇Gx,tr) -γₙ(Gn,tr)]
        
    else
        int = -[nt×(∇G∇B)+κ^2*nt×(Gr) -nt×(∇Gxn) ;
        nt ⋅(∇Gx) -nt ⋅(Gn)]
    end
    return BEAST.matrix_to_bilform(int)
end

Z(κ;rows = [1,2],cols = [1,2],kwargs...) = BEAST.matrix_to_bilform([Bb(κ;kwargs...) T2(κ;kwargs...);T1(κ;kwargs...) A(κ;kwargs...)][rows,cols];kwargs...)
Zp(κ;rows = [1,2],cols = [1,2],kwargs...) = BEAST.matrix_to_bilform([A(κ;kwargs...) T1(κ;kwargs...);T2(κ;kwargs...) Bb(κ;kwargs...)][rows,cols];kwargs...)


function excitation_dirichlet(A,curlA,divA)
    out = [((n × FunctionExcitation{ComplexF64}(A))),
    FunctionExcitation{ComplexF64}(divA)]
    return BEAST.array_to_linform(out)
end
function excitation_neumann(A,curlA,divA)
    out = [((n × FunctionExcitation{ComplexF64}(curlA))),
    BEAST.NDotTrace{ComplexF64}(FunctionExcitation{ComplexF64}(A))]
    return BEAST.array_to_linform(out)
end
function solve(Z,b,X;strat,kwargs...)
    if strat ==:LU
        return BEAST.solve(b,Z,X;kwargs...)
    elseif strat == :GMRES
        return BEAST.gmres_ch(b,Z,X;kwargs...)
    end
    @error "no existing strategy given, use :LU or :GMRES"
end


parent(ind,tree) = tree[ind]
children(ind,tree) = findall(==(ind),tree)


function plot_eigenvalues(M)
    eigs = eigen(M).values
    display(scatter(real.(eigs),imag.(eigs)))

end


### physical constants
ϵ0 = 8.854e-12
μ0 = 4π*1e-7
ω = 10.0^8#*2*pi
κ0 = ω*sqrt(ϵ0*μ0)


#### define meshes
hh = 0.3
Γ1 = meshcuboid(1.0, 1.0, 1.0, hh)
Γ2 =  (@SVector [-1.0,0.0,0.0]) + BEAST.TraceMesh(-Mesh([point(-x,y,z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))))
Γ3 = (@SVector [0.0,0.0,-1.0]) + BEAST.TraceMesh(-translate(Mesh([point(x,y,-z) for (x,y,z) in vertices(Γ1)], deepcopy(cells(Γ1))),[0.0,0.0,0.0]))

Γ = [Γ1,Γ2,Γ3]
Tree = [0,0,0] # give index in which volume material is

HOM = [1,2] #indices of homogeneous domains (without free space)
HomPars = Dict(0=>(ϵ0,μ0),1=>(ϵ0*2,μ0*2),2=>(ϵ0*2,μ0*2))#
κ = Dict(i => ω*sqrt(prod(HomPars[i])) for i in HOM)
κ[0]=κ0
EFIE = [] #indices of pec domains modeled with efie
MFIE = [] #indices of pec domains modeled with mfie
CFIE = [3] #indices of pec domains modeled with cfie
α = Dict(3=>0.2)#index --> alpha

#### Incident field
A(x) = x[1]*[0.0,0.0,1.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0*x[3])
curlA(x) = -[0.0,1.0,0.0]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
divA(x) = -im*κ0 *x[1]*sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])
graddivA(x) = -im*κ0 *sqrt(ϵ0*μ0)*exp(-1.0im*κ0 *x[3])*[1.0,0.0,-im*κ0*x[1]]
#### spaces
Xdb = Γ -> BEAST.DirectProductSpace([raviartthomas(Γ),lagrangec0d1(Γ)])
Xnb = Γ -> BEAST.DirectProductSpace([raviartthomas(Γ),lagrangecxd0(Γ)])
Ydb = Γ -> BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangec0d1(Γ)])
Ynb = Γ -> BEAST.DirectProductSpace([BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ)])

Xdt = Γ -> BEAST.DirectProductSpace([n×raviartthomas(Γ),lagrangecxd0(Γ)])
Xnt = Γ -> BEAST.DirectProductSpace([n×(n×(n×raviartthomas(Γ))),lagrangec0d1(Γ)])
Ydt = Γ -> BEAST.DirectProductSpace([n×BEAST.buffachristiansen(Γ),duallagrangecxd0(Γ)])
Ynt = Γ -> BEAST.DirectProductSpace([n×(n×(n×BEAST.buffachristiansen(Γ))),duallagrangec0d1(Γ)])

Xinn = BEAST.array_to_linform([excitation_neumann(A,curlA,divA)])
Xind = BEAST.array_to_linform([excitation_dirichlet(A,curlA,divA)])
Xin = BEAST.array_to_linform([excitation_neumann(A,curlA,divA),excitation_dirichlet(A,curlA,divA)])

idnd = BEAST.matrix_to_bilform(diagm([Identity(),Identity()]))
@hilbertspace t1 t2
@hilbertspace b1 b2
id = idnd[t1,b1]+idnd[t2,b2]

N = length(Γ)
Q = Dict(i=>diagm([diagm(@SVector [1.0,1.0]),diagm(SVector{2,Float64}([HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]]))]) for i in HOM)
Qp = Dict(i=>diagm([diagm(SVector{2,Float64}([HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]])),diagm(@SVector [1.0,1.0])]) for i in HOM)
Qinv = Dict(i=>diagm([diagm(@SVector [1.0,1.0]),diagm(SVector{2,Float64}(1.0./[HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]]))]) for i in HOM)
Qpinv = Dict(i=>diagm([diagm(SVector{2,Float64}(1.0./[HomPars[parent(i,Tree)][2]/HomPars[i][2],HomPars[i][1]/HomPars[parent(i,Tree)][1]])),diagm(@SVector [1.0,1.0])]) for i in HOM)
t = BEAST.hilbertspace(:t, length(Γ))
b = BEAST.hilbertspace(:b, length(Γ))

##### define space
perm = sortperm([HOM...,EFIE...,CFIE...,MFIE...])
Xt = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnt(Γ[i]),Xdt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xdt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Yt = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ydt(Γ[i]),Ynt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ydt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Xb = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xdb(Γ[i]),Xnb(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnb(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Yb = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynb(Γ[i]),Ydb(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynb(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
Xtmfie = [BEAST.DirectProductSpace[BEAST.DirectProductSpace([Xnt(Γ[i]),Xdt(Γ[i])]) for i in HOM];BEAST.DirectProductSpace[BEAST.DirectProductSpace([Ynt(Γ[i])]) for i in [EFIE...,CFIE...,MFIE...]]][perm]
##### define equation

eqs1 = BEAST.Equation[(Qp[i]*Z(κ[i];tr=-1)*(Qinv[i]))[t[i],b[i]] +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];tr=-1))[t[i],b[j]] for j in HOM ∩ children(i,Tree)]) +
        -sum(BEAST.BilForm[(Qp[i]*Z(κ[i];cols=[2],tr=-1))[t[i],b[j]] for j in [EFIE...,CFIE...,MFIE...] ∩ children(i,Tree)]) ==0 for i in HOM]


eqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
        -sum(BEAST.BilForm[Z(κ0;cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin[t[ci]]
        for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ0;rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idnd[t[ci],b[ci]] ==-Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
            -α[ci]*Xind[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Z(κ0;rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ0;rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idnd()[t[ci],b[ci]] == -(1-α[ci])*Xinn[t[ci]]
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end


eqs3hom = begin BEAST.Equation[(Z(κ[i])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
        -sum(BEAST.BilForm[Z(κ[i];cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
        for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
eqs3efie = begin BEAST.Equation[(Z(κ[i];rows=[2])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
eqs3mfie = begin BEAST.Equation[(Z(κ[i];rows=[1])*(Qinv[i]))[t[ci],b[i]]+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -sum(BEAST.BilForm[Z(κ[i];rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -idnd()[t[ci],b[ci]] ==0
            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ MFIE]end
eqs3cefie = begin BEAST.Equation[α[ci]*(Z(κ[i];rows=[2])*(Qinv[i]))[t[ci],b[i]]+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[2])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[2],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
eqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Z(κ[i];rows=[1])*(Qinv[i]))[t[ci],b[i]]+
            -(1-α[ci])* sum(BEAST.BilForm[Z(κ[i];rows=[1])[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
            -(1-α[ci])*sum(BEAST.BilForm[Z(κ[i];rows=[1],cols=[2])[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) +
            -(1-α[ci])*idnd()[t[ci],b[ci]] == 0
            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

## sum the equations in the two parts in a pmchwt fassion,

symeq = -sum(eqs1)+sum(eqs2cefie)+sum(eqs2efie)+sum(eqs2hom)+sum(eqs3cefie)+sum(eqs3efie)+sum(eqs3hom)
asymeq = sum(eqs2cmfie)+sum(eqs2mfie)+sum(eqs3cmfie)+sum(eqs3mfie)

symfilled = typeof(symeq) <: BEAST.Equation
asymfilled = typeof(asymeq) <: BEAST.Equation


symfilled && (Dsymeq = BEAST.discretise(symeq, (t.∈Xt)..., (b.∈Xb)...))
asymfilled && (Dasymeq = BEAST.discretise(asymeq, (t.∈Xtmfie)..., (b.∈Xb)...))
#assemble system



qs(::BEAST.LocalOperator, a, b) = BEAST.SingleNumQStrat(qsZ[1])
qs(op::BEAST.ComposedOperatorLocal,testspace,trialpsace) = BEAST.SingleNumQStrat(qsZ[2])
qs(op::BEAST.ComposedOperatorIntegral,testspace,trialspace) = BEAST.DoubleNumSauterQstrat(qsZ[3:end]...) 
qs(fn::BEAST.Functional, basis) = BEAST.SingleNumQStrat(8)
asymfilled && ((Za,ba,xx,yy) = assemble(Dasymeq;quadstratfunction = qs))
symfilled && ((Zs,bs,xx,yy) = assemble(Dsymeq;quadstratfunction = qs))

Zunprec = Za*I+Zs*I
#unprec system


bunprec = ba .+ bs 

uz = solve(Zunprec,bunprec,BEAST.DirectProductSpace(Xb);strat=solvestrat,tol=abstol_Z, maxiter=maxiter_Z, restart=restart_Z)
save_uZ && (out[:uZ] = uz[1])
save_uZ && (out[:basis] = Xb)
save_itterZ && (out[:itter_Z] = uz[2])



# ##### nearfield
# using BlockArrays
# function nearfield_B(u,X,κ,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
#     G = BEAST.greenhh3d(wavenumber=κ)
#     gradG = BEAST.∇(G)
#     ∇Gx =  gradG×B
#     Gn = G*(nb*B)
#     #Gnx = G*(n × B)
#     ∇G = gradG*B
#     ∇Gdotn = nb ⋅ (gradG*B)
#     ∇Gdot = B ⋅ gradG
    
#     Gr = G*B
#     ∇G∇B = gradG*BEAST.div(B)
#     ∇Gxn = gradG×(nb*B)
    

#     if HOM
#         @hilbertspace T D R N
#         Bb = -Q[1,1]*potential(∇G∇B, points, u[T], X.factors[1])-Q[1,1]*κ^2* potential(Gr, points, u[T], X.factors[1])
#         Bb .+= Q[2,2]*potential(∇Gxn, points, u[D], X.factors[2])
#         Bb .+= -Q[3,3]*potential(∇Gx, points, u[R], X.factors[3])
        
#     else
#         @hilbertspace R N
#         Bb = -potential(∇Gx, points, u[R], X.factors[1])
        
#     end

#     return sign*Bb
# end

# function nearfield_A(u,X,κ,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
#     G = BEAST.greenhh3d(wavenumber=κ)
#     gradG = BEAST.∇(G)
#     ∇Gx =  gradG×B
#     Gn = G*(nb*B)
#     #Gnx = G*(n × B)
#     ∇G = gradG*B
#     ∇Gdotn = nb ⋅ (gradG*B)
#     ∇Gdot = B ⋅ gradG
    
#     Gr = G*B
#     ∇G∇B = gradG*BEAST.div(B)
#     ∇Gxn = gradG×(nb*B)
    

#     if HOM
#         @hilbertspace T D R N
#         A = -Q[1,1]*potential(∇Gx, points, u[T], X.factors[1])
#         A .+= Q[2,2]*potential(Gn, points, u[D], X.factors[2])
#         A .+= -Q[3,3]*potential(Gr, points, u[R], X.factors[3])
#         A .+= Q[4,4]*potential(∇G, points, u[N], X.factors[4])
#     else
#         @hilbertspace R N
#         A = -potential(Gr, points, u[R], X.factors[1])
#         A .+= potential(∇G, points, u[N], X.factors[2])
#     end

#     return sign*A
# end
# function nearfield_phi(u,X,κ,ω,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
#     G = BEAST.greenhh3d(wavenumber=κ)
#     gradG = BEAST.∇(G)
#     ∇Gx =  gradG×B
#     Gn = G*(nb*B)
#     #Gnx = G*(n × B)
#     ∇G = gradG*B
#     ∇Gdotn = nb ⋅ (gradG*B)
#     ∇Gdot = B ⋅ gradG
    
#     Gr = G*B
#     ∇G∇B = gradG*BEAST.div(B)
#     ∇Gxn = gradG×(nb*B)
    

#     if HOM
#         @hilbertspace T D R N
#         A = Q[2,2]*potential(∇Gdotn, points, u[D], X.factors[2])
#         A .+= -Q[3,3]*potential(∇Gdot, points, u[R], X.factors[3])
#         A .+= -κ^2* Q[4,4]*potential(Gr, points, u[N], X.factors[4])
#     else
#         @hilbertspace R N
#         A = -potential(∇Gdot, points, u[R], X.factors[1])
#         A .+= -κ^2* potential(Gr, points, u[N], X.factors[2])
#     end

#     return sign*A/κ^2*im*ω
# end
# function nearfield_E(u,X,κ,ω,points;sign=1,HOM=true,Q=diagm([1.0,1.0,1.0,1.0]))#nearfield of single volume
#     G = BEAST.greenhh3d(wavenumber=κ)
#     gradG = BEAST.∇(G)
#     graddivG = BEAST.graddiv(G)
#     ∇Gx =  gradG×B
#     Gn = G*(nb*B)
#     #Gnx = G*(n × B)
#     ∇G = gradG*B
#     ∇Gdotn = nb ⋅ (gradG*B)
#     graddivGdotn = graddivG*(nb*B)
#     graddivGdot = graddivG*B
#     ∇Gdot = B ⋅ gradG
    
#     Gr = G*B
#     ∇G∇B = gradG*BEAST.div(B)
#     ∇Gxn = gradG×(nb*B)
    

#     if HOM
#         @hilbertspace T D R N
#         p = Q[2,2]*potential(graddivGdotn, points, u[D], X.factors[2])
#         p .+= -Q[3,3]*potential(graddivGdot, points, u[R], X.factors[3])
#         p .+= -κ^2* Q[4,4]*potential(∇G, points, u[N], X.factors[4])
#     else
#         @hilbertspace R N
#         p = -potential(graddivGdot, points, u[R], X.factors[1])
#         p .+= -κ^2* potential(∇G, points, u[N], X.factors[2])
#     end
#     if HOM
#         @hilbertspace T D R N
#         A = -Q[1,1]*potential(∇Gx, points, u[T], X.factors[1])
#         A .+= Q[2,2]*potential(Gn, points, u[D], X.factors[2])
#         A .+= -Q[3,3]*potential(Gr, points, u[R], X.factors[3])
#         A .+= Q[4,4]*potential(∇G, points, u[N], X.factors[4])
#     else
#         @hilbertspace R N
#         A = -potential(Gr, points, u[R], X.factors[1])
#         A .+= potential(∇G, points, u[N], X.factors[2])
#     end
#     return -sign*p/κ^2*im*ω-sign*im*ω*A
# end
# Xs = range(-2.0,2.0,length=150)
# Zs = range(-1.5,2.5,length=100)
# pts = [point(x,0.5,z) for z in Zs, x in Xs]

# A1 = [nearfield_A(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
# A2 = [nearfield_A(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
# A3 = [nearfield_A(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
# Atot = sum(A1)+sum(A2)+sum(A3)-A.(pts)

# B1 = [nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
# B2 = [nearfield_B(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
# B3 = [nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
# Btot = sum(B1)+sum(B2)+sum(B3)-curlA.(pts)

# H1 = [1/HomPars[parent(i,Tree)][2]*nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts) for i in HOM]
# H2 = [1/HomPars[i][2]*nearfield_B(u[Block(i)],X[i],κ[i],pts;sign=-1,Q=Q[i]^-1) for i in HOM]
# H3 = [1/HomPars[parent(i,Tree)][2]*nearfield_B(u[Block(i)],X[i],κ[parent(i,Tree)],pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
# Htot = sum(H1)+sum(H2)+sum(H3)-1/μ0*curlA.(pts)

# E1 = [nearfield_E(u[Block(i)],X[i],κ[parent(i,Tree)],ω,pts) for i in HOM]
# E2 = [nearfield_E(u[Block(i)],X[i],κ[i],ω,pts;sign=-1,Q=Q[i]^-1) for i in HOM]
# E3 = [nearfield_E(u[Block(i)],X[i],κ[parent(i,Tree)],ω,pts;HOM=false) for i in [EFIE...,MFIE...,CFIE...]]
# Etot = sum(E1)+sum(E2)+sum(E3)-(-im*ω*A.(pts)-graddivA.(pts)/κ0^2*im*ω)

# ##### complement error, different from aps paper, traces from inside are compared, for the reset same formula as in aps paper, for pec all are treated with mfie, as efie yields zero traces in denominator.


# ### to generate plot 
# using Plots
# using LaTeXStrings
# using ColorSchemes
# p1 = Plots.heatmap(Xs, Zs, (real.(getindex.(Htot,2)))*1000,clim=(-0.01001*1000,0.01*1000),show=true,

# xlabel="z-axis [m]",
# ylabel="x-axis [m]",
# colorbar_title=" \nHy [mA/m]",
# c= :PuBu,
# right_margin = 7Plots.mm) 
# rectangle(w, h, x, y) = Plots.Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
# # plot!(rectangle(1,1,0,0),fillcolor=plot_color(nothing),linecolor=:black,linewidth=4)
# # plot!(rectangle(1,1,-1,0),fillcolor=plot_color(nothing),linecolor=:black,linewidth=4)
# # plot!(rectangle(1,1,0,-1),fillcolor=plot_color(nothing),linecolor=:black,linewidth=4,legend=false)
# # annotate!(-0.5,0.5,"Ω₃", :black)
# # annotate!(0.5,0.5,"Ω₁", :black)
# # annotate!(0.5,-0.5,"Ω₂", :black)
# display(plot!(
#     xtickfont=font(11),
#     ytickfont=font(11),
#     xguidefontsize=12,
# yguidefontsize=13,
# colorbar_titlefontsize=13
# ))
# #savefig("Hfield.pdf")



# ##### define equation
# ceqs2hom = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#            -sum(BEAST.BilForm[Z(κ0;cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xin_notdual[t[ci]]
#            for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
# # ceqs2efie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
# #             -sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xind_notdual[t[ci]]
# #             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
# ceqs2mfie = begin BEAST.Equation[-sum(BEAST.BilForm[Z(κ0;rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==-Xinn_notdual[t[ci]]
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ [MFIE...,EFIE...]]end
# # ceqs2cefie =begin  BEAST.Equation[-α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
# #             -α[ci]*sum(BEAST.BilForm[Z(κ0;rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==
# #             -α[ci]*Xind_notdual[t[ci]] for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
# ceqs2cmfie = begin BEAST.Equation[-(1-α[ci])* sum(BEAST.BilForm[Z(κ0;rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -(1-α[ci])*sum(BEAST.BilForm[Z(κ0;rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)])  == -(1-α[ci])*Xinn_notdual[t[ci]]
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end

# # deqs2hom = BEAST.discretise(eqs2hom, (t.∈X),(b.∈X))    
# # deqs2efie = BEAST.discretise(eqs2efie, (t.∈X),(b.∈X))    
# # deqs2mfie = BEAST.discretise(eqs2mfie, (t.∈Y),(b.∈X))    
# # deqs2cefie = BEAST.discretise(eqs2cefie, (t.∈X),(b.∈X))    
# # deqs2cmfie = BEAST.discretise(eqs2cmfie, (t.∈Y),(b.∈X))    


# ##### define equation
# ceqs3hom = begin BEAST.Equation[(Z(κ[i];tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
#             -sum(BEAST.BilForm[Z(κ[i];tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#            -sum(BEAST.BilForm[Z(κ[i];cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
#            for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ HOM] end
# # ceqs3efie = begin BEAST.Equation[(Z(κ[i];rows=[3,4],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
# #             -sum(BEAST.BilForm[Z(κ[i];rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
# #             -sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
# #             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ EFIE] end
# ceqs3mfie = begin BEAST.Equation[(Z(κ[i];rows=[1,2],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
#             -sum(BEAST.BilForm[Z(κ[i];rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)])  ==0
#             for i in [0], ci in 1:N if ci ∈ children(i,Tree) ∩ [MFIE...,EFIE...]]end
# # ceqs3cefie = begin BEAST.Equation[α[ci]*(Z(κ[i];rows=[3,4],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
# #             -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
# #             -α[ci]*sum(BEAST.BilForm[Z(κ[i];rows=[3,4],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) ==0
# #             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end
# ceqs3cmfie = begin BEAST.Equation[(1-α[ci])*(Z(κ[i];rows=[1,2],tr=-1,dual=false)*(Q[i]^-1))[t[ci],b[i]]+
#             -(1-α[ci])* sum(BEAST.BilForm[Z(κ[i];rows=[1,2],tr=-1,dual=false)[t[ci],b[j]] for j in HOM ∩ children(i,Tree)])+
#             -(1-α[ci])*sum(BEAST.BilForm[Z(κ[i];rows=[1,2],cols=[3,4],tr=-1,dual=false)[t[ci],b[j]] for j in [EFIE...,MFIE...,CFIE...] ∩ children(i,Tree)]) == 0
#             for i in HOM, ci in 1:N if ci ∈ children(i,Tree) ∩ CFIE] end



# ceq = sum(ceqs2hom)+sum(ceqs3hom)
# cmeq = sum(ceqs2cmfie)+sum(ceqs2mfie)+sum(ceqs3cmfie)+sum(ceqs3mfie)

# Dseq = BEAST.discretise(ceq+cmeq, (t.∈X)..., (b.∈X)...)

# F,bsa,_,_ = assemble(Dseq)

# using NestedUnitRanges

# ax = nestedrange(BEAST.DirectProductSpace(X), 1,numfunctions)
# u_error = PseudoBlockVector(F*u-bsa, (ax,))
# #select in u_error en in u

# function select_trace(u,trace_ind)
#     usub = 0*deepcopy(u)
#     N = length(blocks(u))
#     println(N)
#     for i in 1:N
#         if i ∈ HOM
            
#             usub2 = deepcopy(usub[Block(i)])
#             usub2[Block(trace_ind)] .+= u[Block(i)][Block(trace_ind)] 
#             usub[Block(i)] .+= usub2
#         elseif i∈ [EFIE...,MFIE...,CFIE...] && trace_ind ∈ [3,4]
#             println("pec")
#             usub2 = deepcopy(usub[Block(i)])
#             usub2[Block(trace_ind-2)] .+= u[Block(i)][Block(trace_ind-2)]
#             usub[Block(i)] .+= usub2
#         end
#     end
#     return usub
# end

# Gsym = assemble(BEAST.diag(BEAST.diag(Identity())),BEAST.DirectProductSpace(X),BEAST.DirectProductSpace(X))
# Ginvsym = inv(Matrix(Gsym))

# nom = sqrt.([Array(select_trace(u_error,i))'*Ginvsym*Array(select_trace(u_error,i))  for i in 1:4])
# denom = sqrt.([Array(select_trace(u,i))'*Matrix(Gsym)*Array(select_trace(u,i))  for i in 1:4])
# error = 1/4*sum(nom./denom)* 