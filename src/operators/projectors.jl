abstract type QHComponent end

#Local loops + global loops
struct Loops <: QHComponent end
#Stars
struct Stars <: QHComponent end
#Stars
struct DualStars <: QHComponent end
#Local loops + global loops
struct DualLoops <: QHComponent end

abstract type ComputeStrat end
struct Direct <: ComputeStrat end
struct IterativePackage{T} <: ComputeStrat end
const Iterative = IterativePackage{:IterativeSolvers}

struct QHProjector{C<:QHComponent, S<:ComputeStrat} 
    kwargs::Dict{Symbol,Any}
end

function PΣ(;compStrat = Iterative, kwargs...)
    return QHProjector{Stars,compStrat}(kwargs)
end

function PΛ(;compStrat = Iterative, kwargs...)
    return QHProjector{Loops,compStrat}(kwargs)
end

function ℙΣ(;compStrat = Iterative, kwargs...)
    return QHProjector{DualLoops,compStrat}(kwargs)
end

function ℙΛ(;compStrat = Iterative, kwargs...)
    return QHProjector{DualStars,compStrat}(kwargs)
end

"""
saddlepoint(A,B,P1,P2)
Create preconditioned saddlepoint matrix from A and B with block diagonal preconditoners [P1,P2]1 and P2. 
"""
function saddlepoint(method,A::SparseMatrixCSC,B::SparseMatrixCSC,P1::SparseMatrixCSC,P2::SparseMatrixCSC; kwargs...) 
    T= eltype(A)
    nP = size(B,2)
    SP = [A    B
          B'  spzeros(T,nP,nP)]
    Pdiv = blockdiag(CholeskyFactorization(P1),CholeskyFactorization(P2))
    return create_it_projector(method, SP;left_preconditioner=Pdiv, kwargs...)
end
function create_it_projector(::Type{IterativePackage{:IterativeSolvers}}, SP;left_preconditioner= nothing, maxiter = 200, restart = 50,
    reltol = 1e-8, verbose = false, kwargs...)
    return GMRESSolver(SP, left_preconditioner=left_preconditioner, maxiter=maxiter, restart=restart, reltol=reltol, verbose=verbose; kwargs...)
end
function create_it_projector(::Type{IterativePackage{:Krylov}}, SP;left_preconditioner = LinearAlgebra.I, restart=true, itmax = 200, type=eltype(SP), 
    memory = 50, workspace=GMRESWorkspace(SP; type=type, memory=memory), verbose = false, kwargs...)
    M = left_preconditioner
    return GMRES(SP;M=M, restart=restart, itmax=itmax,memory = memory, workspace=workspace, kwargs...)
end
function assemble(::QHProjector, X::Space; quadstrat=defaultquadstrat)
    error("Not implemented")
end

#RT Basis
function assemble(::QHProjector{Stars,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    Σ = Matrix(connectivity(X.geo, edges, sign))
    return Σ * pinv(Σ'*Σ) * Σ'
end

function assemble(::QHProjector{Loops,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
   edges = setminus(skeleton(X.geo,1), boundary(X.geo))
   Σ = Matrix(connectivity(X.geo, edges, sign))
   return LinearAlgebra.I - Σ*pinv(Σ'*Σ)*Σ'
end

function assemble(::QHProjector{DualStars,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    verts = setminus(skeleton(X.geo,0), skeleton(boundary(X.geo),0))
    Λ = Matrix(connectivity(verts, edges, sign))
    return  Λ * pinv(Λ'*Λ) * Λ'
end

function assemble(::QHProjector{DualLoops,Direct}, X::RTBasis; quadstrat=defaultquadstrat)
    edges = setminus(skeleton(X.geo,1), boundary(X.geo))
    verts = setminus(skeleton(X.geo,0), skeleton(boundary(X.geo),0))
    Λ = Matrix(connectivity(verts, edges, sign))
    return  LinearAlgebra.I - Λ * pinv(Λ'*Λ) * Λ'
end

function assemble(QHP::QHProjector{Stars,IterativeMethod}, X::RTBasis; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    #create auxilarry basis functions
    P = lagrangecxd0(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilarry matrices
    # Gxx = assemble(Identity(),X,X;quadstrat)
    T = scalartype(X)
    Gxx = sparse(T,I,nX,nX)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,Gpp;QHP.kwargs...)
    return P0Σ*SP*Px0
end

function assemble(::QHProjector{Loops,IterativeMethod}, X::RTBasis; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    #create auxilarry basis functions
    P = lagrangecxd0(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilary matrices
    # Gxx = assemble(Identity(),X,X;quadstrat)
    T = scalartype(X)
    Gxx = sparse(T,I,nX,nX)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,Gpp)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualLoops,IterativeMethod}, X::RTBasis; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    #create auxilarry basis functions
    P = duallagrangecxd0(X.geo)
    X = buffachristiansen(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)

    T = scalartype(X)
    Gxx = sparse(T,I,nX,nX)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    sqGpp = sqrt.(Gpp)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)*sqGpp
    Px0 = [Gxx
           spzeros(nP,nX)]
    # P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,sqGpp*Gpp*sqGpp)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualStars,IterativeMethod}, X::RTBasis; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}

    P = duallagrangecxd0(X.geo)
    X = buffachristiansen(X.geo)
    nX = numfunctions(X)
    nP = numfunctions(P)

    T = scalartype(X)
    Gxx = sparse(T,I,nX,nX)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    sqGpp = sqrt.(Gpp)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)*sqGpp
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,sqGpp*Gpp*sqGpp)
    return P0Σ*SP*Px0
end
#GWP basis
function assemble(::QHProjector{Stars,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    return Σp * pinv(Matrix(Σp'*inv(Matrix(Gxx))*Σp)) * Σp'
end

function assemble(::QHProjector{Loops,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    return Gxx - Σp * pinv(Matrix(Σp'*inv(Matrix(Gxx))*Σp)) * Σp'
end

function assemble(::QHProjector{DualLoops,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat) 
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    return Gxx - Λp * pinv(Matrix(Λp'*inv(Matrix(Gxx))*Λp)) * Λp'
end

function assemble(::QHProjector{DualStars,Direct}, X::GWPDivSpace; quadstrat=defaultquadstrat)
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat) 
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    return Λp * pinv(Matrix(Λp'*inv(Matrix(Gxx))*Λp)) * Λp'
end

function assemble(::QHProjector{Stars,IterativeMethod}, X::GWPDivSpace; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    p = X.degree
    #create auxilarry basis functions
    P = lagrangecx(X.geo,order=p)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Σ = [spzeros(nX,nX) Σp]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,Gpp)
    return P0Σ*SP*Px0
end

function assemble(::QHProjector{Loops,IterativeMethod}, X::GWPDivSpace; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    p = X.degree
    #create auxilary basis functions
    P = lagrangecx(X.geo,order=p)
    nX = numfunctions(X)
    nP = numfunctions(P)
    #assemble auxilary matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Dxx = assemble(Identity(),divergence(X),divergence(X);quadstrat)
    Gpp = assemble(Identity(),P,P;quadstrat)
    Σp = assemble(Identity(),divergence(X),P;quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    SP = saddlepoint(IterativeMethod, Gxx,Σp,Gxx+Dxx,Gpp)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualLoops,IterativeMethod}, X::GWPDivSpace; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    nX = numfunctions(X)
    nL = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    #Σp = assemble(Identity(),divergence(X),P)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nL,nX)]
    SP = saddlepoint(IterativeMethod, Gxx,Λp,Gxx,Gll+Cll)
    return Px0'*SP*Px0
end

function assemble(::QHProjector{DualStars,IterativeMethod}, X::GWPDivSpace; quadstrat=defaultquadstrat) where {IterativeMethod <:IterativePackage}
    p = X.degree
    #create auxilarry basis functions
    L = lagrangec0(X.geo,order=p+1)
    nX = numfunctions(X)
    nP = numfunctions(L)
    #assemble auxilarry matrices
    Gxx = assemble(Identity(),X,X;quadstrat)
    Cll = assemble(Identity(),curl(L),curl(L);quadstrat)
    Gll = assemble(Identity(),L,L;quadstrat)
    #Σp = assemble(Identity(),divergence(X),P)
    Λp = assemble(Identity(),X,curl(L);quadstrat)
    Px0 = [Gxx
           spzeros(nP,nX)]
    P0Λ = [spzeros(nX,nX) Λp]
    SP = saddlepoint(IterativeMethod, Gxx,Λp,Gxx,Gll+Cll)
    return P0Λ*SP*Px0
end   


@testitem "QH-Projectors" begin 

    using CompScienceMeshes, BEAST
    using LinearAlgebra

    h = 0.5
    M = meshrectangle(1.0,1.0,h)

    X = raviartthomas(M)
    Y = BEAST.gwpdiv(M;order=1)

    G = assemble(BEAST.Identity(),Y,Y)
    iG = BEAST.cholesky(G)

    b = rand(numfunctions(Y))
    d = rand(numfunctions(X))

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Direct),X)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Direct),X)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Direct),X)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Direct),X)

    @test norm(PΣd*PΛd*d)/norm(d) < sqrt(eps())
    @test norm(ℙΣd*ℙΛd*d)/norm(d) < sqrt(eps())
≈
    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Iterative),X)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Iterative),X)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Iterative),X)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Iterative),X)

    @test norm(PΣd*PΛd*d)/norm(d) < sqrt(eps())
    @test norm(ℙΣd*ℙΛd*d)/norm(d) < sqrt(eps())

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Direct),Y)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Direct),Y)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Direct),Y)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Direct),Y)

    @test norm(PΣd*iG*PΛd*b)/norm(b) <  sqrt(eps())
    @test norm(ℙΣd*iG*ℙΛd*b)/norm(b) <  sqrt(eps())

    PΣd = assemble(BEAST.PΣ(;compStrat =BEAST.Iterative),Y)
    PΛd = assemble(BEAST.PΛ(;compStrat =BEAST.Iterative),Y)

    ℙΣd = assemble(BEAST.ℙΣ(;compStrat =BEAST.Iterative),Y)
    ℙΛd = assemble(BEAST.ℙΛ(;compStrat =BEAST.Iterative),Y)

    @test norm(PΣd*iG*PΛd*b)/norm(b) <  sqrt(eps())
    @test norm(ℙΣd*iG*ℙΛd*b)/norm(b) <  sqrt(eps())

end 