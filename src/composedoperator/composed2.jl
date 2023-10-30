import Base: *, div
import LinearAlgebra: ×, ⋅




#const i4pi = 1 / (4pi)
abstract type ComposedOperatorLocal <: LocalOperator end
abstract type ComposedOperatorIntegral <: IntegralOperator end
abstract type Kernel <: ComposedOperatorIntegral end
ComposedOperator = Union{ComposedOperatorIntegral,ComposedOperatorLocal}

struct BasisFunction <: ComposedOperatorLocal end
struct DivBasisFunction <: ComposedOperatorLocal end
const B = BasisFunction
export B


struct TestNormal <: ComposedOperatorLocal end
struct TrialNormal <: ComposedOperatorLocal end
struct TraceDirection <: ComposedOperatorLocal end
scalartype(op::Union{TestNormal,TrialNormal,BasisFunction,TraceDirection}) = Float16


const nt = TestNormal()

struct Potential{T} <: AbstractOperator
    operator::T
    surface
end
normal_surface(p::Potential) = p.surface
function Potential(linop::LinearCombinationOfOperators,surface)
    out = ZeroOperator()
    for (op,coeff) in zip(linop.ops,linop.coeffs)
        out += coeff*Potential(op,surface)
    end
    return out
end
struct TraceOperator{T} <: AbstractOperator
    operator::T
    directionsurface
    normalsurface
end
function TraceOperator(linop::LinearCombinationOfOperators,dsurf,nsurf)
    out = ZeroOperator()
    for (op,coeff) in zip(linop.ops,linop.coeffs)
        out += coeff*TraceOperator(op,dsurf,nsurf)
    end
return out
end
normal_surface(p::TraceOperator) = p.normalsurface
direction_surface(p::TraceOperator) = p.directionsurface
struct TimesLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct DotLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct CrossLocal{U,V} <: ComposedOperatorLocal
    lhs::U
    rhs::V
end

struct TimesIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

struct DotIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

struct CrossIntegral{U,V} <: ComposedOperatorIntegral
    lhs::U
    rhs::V
end

function TimesLocal(lhs::LinearCombinationOfOperators,rhs::Union{ComposedOperator,LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*TimesLocal(op,rhs)
    end
    return out
end

function TimesLocal(lhs::ComposedOperator,rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*TimesLocal(lhs,op)
    end
    return out
end

function DotLocal(lhs::LinearCombinationOfOperators,rhs::Union{ComposedOperator,LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*DotLocal(op,rhs)
    end
    return out
end

function DotLocal(lhs::ComposedOperator,rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*DotLocal(lhs,op)
    end
    return out
end

function CrossLocal(lhs::LinearCombinationOfOperators,rhs::Union{ComposedOperator,LinearCombinationOfOperators})
    out = ZeroOperator()
    for (op,coeff) in zip(lhs.ops,lhs.coeffs)
        out += coeff*CrossLocal(op,rhs)
    end
    return out
end

function CrossLocal(lhs::ComposedOperator,rhs::LinearCombinationOfOperators)
    out = ZeroOperator()
    for (op,coeff) in zip(rhs.ops,rhs.coeffs)
        out += coeff*CrossLocal(lhs,op)
    end
    return out
end
get_constructor(::TimesLocal) = TimesLocal
get_constructor(::TimesIntegral) = TimesLocal
get_constructor(::DotLocal) = DotLocal
get_constructor(::DotIntegral) = DotLocal
get_constructor(::CrossLocal) = CrossLocal
get_constructor(::CrossIntegral) = CrossLocal

OperationsLocal = Union{TimesLocal,DotLocal,CrossLocal}
OperationsIntegral = Union{TimesIntegral,DotIntegral,CrossIntegral}
Operations = Union{OperationsIntegral,OperationsLocal}

TimesLocal(a::ComposedOperatorIntegral,b::ComposedOperatorLocal) = TimesIntegral(a,b)
TimesLocal(a::ComposedOperatorLocal,b::ComposedOperatorIntegral) = TimesIntegral(a,b)
TimesLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = TimesIntegral(a,b)

DotLocal(a::ComposedOperatorIntegral,b::ComposedOperatorLocal) = DotIntegral(a,b)
DotLocal(a::ComposedOperatorLocal,b::ComposedOperatorIntegral) = DotIntegral(a,b)
DotLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = DotIntegral(a,b)

CrossLocal(a::ComposedOperatorIntegral,b::ComposedOperatorLocal) = CrossIntegral(a,b)
CrossLocal(a::ComposedOperatorLocal,b::ComposedOperatorIntegral) = CrossIntegral(a,b)
CrossLocal(a::ComposedOperatorIntegral,b::ComposedOperatorIntegral) = CrossIntegral(a,b)

×(a::ComposedOperator,b::ComposedOperator) = CrossLocal(a,b)
⋅(a::ComposedOperator,b::ComposedOperator) = DotLocal(a,b)
*(a::ComposedOperator,b::ComposedOperator) = TimesLocal(a,b)

Base.div(::BasisFunction) = DivBasisFunction()

scalartype(op::Operations) = promote_type(scalartype(op.lhs),scalartype(op.rhs))
scalartype(op::Potential) = scalartype(op.operator)
scalartype(op::TraceOperator) = scalartype(op.operator)

function count_test_normals(op::Operations)
    count_test_normals(op.lhs) + count_test_normals(op.rhs)
end
count_test_normals(op::Kernel) = 0
count_test_normals(op::TestNormal) = 1
count_test_normals(op::TrialNormal) = 0

function count_trial_normals(op::Operations)
    count_trial_normals(op.lhs) + count_trial_normals(op.rhs)
end
count_trial_normals(op::Kernel) = 0
count_trial_normals(op::TestNormal) = 0
count_trial_normals(op::TrialNormal) = 1

function replace_normal_by_testnormal(op::Operations)
    get_constructor(op)(replace_normal_by_testnormal(op.lhs),replace_normal_by_testnormal(op.rhs))
end
replace_normal_by_testnormal(op::Kernel) = op
replace_normal_by_testnormal(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction}) = op
replace_normal_by_testnormal(op::NormalVector) = TestNormal()

function replace_normal_by_trialnormal(op::Operations)
    get_constructor(op)(replace_normal_by_trialnormal(op.lhs),replace_normal_by_trialnormal(op.rhs))
end
replace_normal_by_trialnormal(op::Kernel) = op
replace_normal_by_trialnormal(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction}) = op
replace_normal_by_trialnormal(op::NormalVector) = TrialNormal()



function build_potential(op::ComposedOperator,surface::CompScienceMeshes.Mesh)
    newop = replace_normal_by_trialnormal(op)
    Potential(newop,surface)
end

function γ(op::Operations,sign)

    return get_constructor(op)(γ(op.lhs,sign),γ(op.rhs,sign))
end
γ(op::Union{TestNormal,TrialNormal,BasisFunction,DivBasisFunction},sign) = op
function check_if_coincide(a,b) 
    @warn "all meshes coincide"
    return true
end

function γ(op::Potential,direction::CompScienceMeshes.AbstractMesh,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
    check_if_coincide(op.surface,surface) || return op.operator
    newop = γ(op.operator,sign)
    return TraceOperator(Potential(newop,op.surface),direction,surface)
end
function γ(op::Potential,surface::CompScienceMeshes.AbstractMesh,sign::Int)# sign + if according to normal on surface, - otherwise
    γ(op,surface,surface,sign)
end


# γₜᶜ(op::Potential,surface) = nt×(γ(op,surface,-1)×nt)
# γₜ(op::Potential,surface) = nt×(γ(op,surface,1)×nt)

# γₛ(op::Potential,surface) = nt×γ(op,surface,1)
# γₛᶜ(op::Potential,surface) = nt×γ(op,surface,-1)
# γₙ(op::Potential,surface) = nt⋅γ(op,surface,1)
# γₙᶜ(op::Potential,surface) = nt⋅γ(op,surface,-1)
# τ(op::Potential,surface) = γ(op,surface,1)
# τᶜ(op::Potential,surface) = γ(op,surface,-1)

function greenhh3d(;
    gamma=nothing,
    wavenumber=nothing)

gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
@assert gamma !== nothing
GreenHH3D(gamma)
end
function gradgreenhh3d(;
    gamma=nothing,
    wavenumber=nothing)

gamma, wavenumber = gamma_wavenumber_handler(gamma, wavenumber)
@assert gamma !== nothing
GradGreenHH3D(gamma)
end

struct GreenHH3D{T} <: Kernel
    gamma::T
end

struct GradGreenHH3D{T} <: Kernel
    gamma::T
end


function (op::GreenHH3D)(x,y,g)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(i4pi*iR)
    Ref(green)
end

function (op::GradGreenHH3D)(x,y,g)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-gamma*R)*(iR*i4pi)
    gradgreen = -(gamma + iR) * green * (iR * r)
    Ref(gradgreen)
end

γ(op::GreenHH3D,sign::Int) = op
γ(op::GradGreenHH3D,sign::Int) = op + sign*1/2*TraceDirection()

grad(G::GreenHH3D) = GradGreenHH3D(G.gamma)

function (::Nabla)(G::Kernel)
    grad(G)
end


scalartype(G::GreenHH3D{T}) where {T} = T
scalartype(G::GradGreenHH3D{T}) where {T} = T


function (op::Union{TimesIntegral,TimesLocal})(x,y,g)
    op.lhs(x,y,g).*op.rhs(x,y,g)
end
function (op::Union{DotIntegral,DotLocal})(x,y,g)
    dot.(op.lhs(x,y,g),op.rhs(x,y,g))
end
function (op::Union{CrossIntegral,CrossLocal})(x,y,g)
    cross.(op.lhs(x,y,g),op.rhs(x,y,g))
end

function (op::TrialNormal)(x,y,g)
    Ref(normal(y))
end
function (op::TestNormal)(x,y,g)
    Ref(normal(x))
end
function (op::TraceDirection)(x,y,g)
    Ref(tracedirection(x))
end
function (op::BasisFunction)(x,y,g)
    getvalue(g)
end
function (op::DivBasisFunction)(x,y,g)
    getdivergence(g)
end



function (igd::Integrand{<:ComposedOperatorIntegral})(x,y,f,g)
    _krondot(getvalue(f),igd.operator(x,y,g))
end

function integrand(op::ComposedOperatorLocal,kernel,x,y,f,g)
    _krondot(getvalue(f),op(x,y,g))
end


defaultquadstrat(op::ComposedOperatorIntegral,testspace::Space,trialspace::Space) = DoubleNumSauterQstrat(6,7,5,5,4,3) 
defaultquadstrat(op::ComposedOperatorLocal,testspace::Space,trialpsace::Space) = SingleNumQStrat(6)
sign_upon_permutation(op::ComposedOperator,I,J) = Combinatorics.levicivita(I)^count_test_normals(op)*Combinatorics.levicivita(J)^count_trial_normals(op)


normalorient(op::ComposedOperator,sign_test_normal,sign_trial_normal) = sign_test_normal^count_test_normals*sign_trial_normal^count_trial_normals


defaultquadstrat(op::TraceOperator,test,trial) = defaultquadstrat(op.operator,test,trial)
defaultquadstrat(op::Potential,test,trial) = defaultquadstrat(op.operator,test,trial)



function assemble!(op::Potential, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi}; 
    quadstrat = defaultquadstrat(op, test_functions, trial_functions))

    nsurf = normal_surface(op)
    surf = geometry(trial_functions)
    @assert same_geometry(nsurf,surf)
    trial = typeof(trial_functions).name.wrapper
    trial_functions = trial(OrientedMesh(surf,nsurf),trial_functions.fns,trial_functions.pos)

    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end

function assemble!(op::TraceOperator, test_functions::Space, trial_functions::Space,
    store, threading = Threading{:multi}; 
    quadstrat = defaultquadstrat(op, test_functions, trial_functions))

    nsurf = normal_surface(op)
    dsurf = direction_surface(op)
    surf = geometry(test_functions)
    @assert same_geometry(nsurf,dsurf)
    @assert same_geometry(nsurf,surf)

    test = typeof(test_functions).name.wrapper

    test_functions = test(TraceMesh(OrientedMesh(surf,nsurf),dsurf),test_functions.fns,test_functions.pos)

    assemble!(op.operator, test_functions, trial_functions, store, threading;
    quadstrat = quadstrat)
end

kernelvals(op::ComposedOperator,a) = nothing