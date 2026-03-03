import Krylov
import Suppressor

struct GMRES{T} <: LinearMap{T}
    A
    M
    N
    ldiv::Bool
    restart::Bool
    reorthogonalization::Bool
    atol::T
    rtol::T
    itmax::Int
    timemax::Float64
    verbose::Int
    history::Bool
    callback::Function
    iostream::IO
    memory::Int
    workspace
end

Base.axes(solver::GMRES) = reverse(axes(solver.A))
Base.size(solver::GMRES) = reverse(size(solver.A))

LinearAlgebra.adjoint(solver::GMRES) = GMRES(adjoint(solver.A); M=solver.N, N=solver.M, memory=solver.memory, restart=solver.restart,
                                              atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax, verbose=solver.verbose)
LinearAlgebra.transpose(solver::GMRES) = GMRES(transpose(solver.A); M=solver.N, N=solver.M, memory=solver.memory, restart=solver.restart,
                                              atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax, verbose=solver.verbose)

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::GMRES, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    solver.verbose > 0 && println("Number of iterations: ", ch.niter)
    ch.solved || error("Iterative solver did not converge.")
    return y
end


function GMRES(A;
    M = LinearAlgebra.I,
    N = LinearAlgebra.I,
    ldiv = false,
    restart = false,
    reorthogonalization = false,
    atol = zero(real(eltype(A))),
    rtol = sqrt(eps(real(eltype(A)))),
    itmax = 2 * size(A,2),
    timemax = Inf,
    verbose = 1,
    history = true,
    callback = x -> false,
    iostream = Krylov.kstdout,
    memory = 50,
    workspace = GMRESWorkspace(A; memory = memory, type = eltype(A)))

    # m, n = size(A)
    T = eltype(A)
    # ws = Krylov.GmresWorkspace(m, n, Vector{T}; memory)
    return GMRES{T}(A, M, N, ldiv, restart, reorthogonalization, atol, rtol, itmax, timemax, verbose, history, callback, iostream, memory,workspace)
end
function GMRESWorkspace(A; memory = 50, type = eltype(A)) 
    m, n = size(A)
    T = type
    ws = Krylov.GmresWorkspace(m, n, Vector{T}; memory)
    return ws
end

GMRESWorkspace(N::Int, type::Type{<:Number}, memory::Int) = Krylov.GmresWorkspace(N, N, Vector{type}; memory)
function solve!(x, solver::GMRES{<: Complex}, b::AbstractVector{<:Complex})

    # m, n = size(solver.A)
    T = promote_type(eltype(solver.A), eltype(b))
    # @show T
    # ws = Krylov.GmresWorkspace(m, n, Vector{T}; solver.memory)
    solver.workspace === nothing ? ws = GMRESWorkspace(solver.A; type = T, memory = solver.memory) : ws = solver.workspace
    # @show solver.workspace

     
        Krylov.gmres!(ws, solver.A, Vector(b);
            M = solver.M, 
            N = solver.N,
            # memory = solver.memory, 
            restart = solver.restart, 
            atol = solver.atol, 
            rtol = solver.rtol, 
            itmax = solver.itmax, 
            verbose = solver.verbose,
            history = solver.history,
            callback = solver.callback,
            iostream = solver.iostream,
            reorthogonalization = solver.reorthogonalization,
            timemax = solver.timemax,
            ldiv = solver.ldiv)
      
    x .= Krylov.solution(ws)

    return x, Krylov.statistics(ws), ws
end
function solve!(x, solver::GMRES{<: Real}, b::AbstractVector{<:Real})

    # m, n = size(solver.A)
    T = promote_type(eltype(solver.A), eltype(b))
    # @show T
    # ws = Krylov.GmresWorkspace(m, n, Vector{T}; solver.memory)
    solver.workspace === nothing ? ws = GMRESWorkspace(solver.A; type = T, memory = solver.memory) : ws = solver.workspace
    # @show solver.workspace

     
        Krylov.gmres!(ws, solver.A, Vector(b);
            M = solver.M, 
            N = solver.N,
            # memory = solver.memory, 
            restart = solver.restart, 
            atol = solver.atol, 
            rtol = solver.rtol, 
            itmax = solver.itmax, 
            verbose = solver.verbose,
            history = solver.history,
            callback = solver.callback,
            iostream = solver.iostream,
            reorthogonalization = solver.reorthogonalization,
            timemax = solver.timemax,
            ldiv = solver.ldiv)
      
    x .= Krylov.solution(ws)

    return x, Krylov.statistics(ws), ws
end

function solve!(x, solver::GMRES{<:Real}, b::AbstractVector{<:Complex})

    # m, n = size(solver.A)
    T = eltype(solver.A)
    # @show T
    # ws = Krylov.GmresWorkspace(m, n, Vector{T}; solver.memory)
    solver.workspace === nothing ? ws = GMRESWorkspace(solver.A; type = T, memory = solver.memory) : ws = solver.workspace
    # @show solver.workspace

    #  Suppressor.@suppress_err begin
        Krylov.gmres!(ws, solver.A, real.(Vector(b));
            M = solver.M, 
            N = solver.N,
            # memory = solver.memory, 
            restart = solver.restart, 
            atol = solver.atol, 
            rtol = solver.rtol, 
            itmax = solver.itmax, 
            verbose = solver.verbose,
            history = solver.history,
            callback = solver.callback,
            iostream = solver.iostream,
            reorthogonalization = solver.reorthogonalization,
            timemax = solver.timemax,
            ldiv = solver.ldiv)
    #   end
    x .= Krylov.solution(ws)
    #  Suppressor.@suppress_err begin
        Krylov.gmres!(ws, solver.A, imag.(Vector(b));
            M = solver.M, 
            N = solver.N,
            # memory = solver.memory, 
            restart = solver.restart, 
            atol = solver.atol, 
            rtol = solver.rtol, 
            itmax = solver.itmax, 
            verbose = solver.verbose,
            history = solver.history,
            callback = solver.callback,
            iostream = solver.iostream,
            reorthogonalization = solver.reorthogonalization,
            timemax = solver.timemax,
            ldiv = solver.ldiv)
    #   end
        x .+= im .* Krylov.solution(ws)
    return x, Krylov.statistics(ws), ws
end

function solve(solver::GMRES, b)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Array{T}, axes(solver)[2])
    fill!(x,0)
    x, ch, ws = solve!(x, solver, b)
    # z = similar(Array{T}, axes(solver)[2])
    # copyto!(z,x)
    return x, ch, ws
end

@testitem "GMRES" begin
    using LinearAlgebra

    A = rand(10, 10)
    b = rand(10)

    Ai = BEAST.GMRES(A; restart=true, atol=1e-8, rtol=1e-8, verbose=1, memory=200)
    x, ch = solve(Ai, b)

    Ai2 = BEAST.GMRESSolver(A; restart=200, abstol=1e-8, reltol=1e-8, verbose=true)
    x2, ch2 = solve(Ai2, b)

    @test norm(x - x2) < 1e-6
end