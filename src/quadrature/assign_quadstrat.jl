struct AssignQuadStrat <: AbstractOperator
    op::AbstractOperator
    qs
end
# struct DefaultQuadstrat end
# AssignQuadStrat(op,qs::DefaultQuadstrat) = op 
# function AssignQuadStrat(op::AssignQuadStrat,qs) 
#     @warn "tryed but failed to overwrite quadstrat "*string(typeof(op.qs))*" assigned to "*string(typeof(op))*" with "*string(typeof(qs))
#     op
# end
# assign_quadstrat(op::Operator, qs::AbstractQuadStrat) = AssignQuadStrat(op, qs)
# operatortype(a::AbstractQuadStrat) = error("assign operatortype for quadstrat "*string(typeof(a)))

# assign_quadstrat(op::LinearCombinationOfOperators, qs::AbstractQuadStrat) = sum([c*assign_quadstrat(opi, qs) for (c,opi) in zip(op.coeffs,op.ops)])
# function assign_quadstrat(op::LinearCombinationOfOperators, qss::Vector{<:AbstractQuadStrat}) 
#     assign_quadstrat(op,qss,operatortype.(qss))
# end
# function assign_quadstrat(op::LinearCombinationOfOperators, qss::Vector{<:AbstractQuadStrat}, operator_types::Vector)
#     return sum([c*assign_quadstrat(opi,qss[findfirst(x-> typeof(opi) ∈ x ,operator_types)]) for (c,opi) in zip(op.coeffs,op.ops)])
# end

# assign_quadstrat(op::BilForm, qs) = BilForm(op.test_space,op.trial_space,assign_quadstrat(op.terms,qs))
# assign_quadstrat(op::BilTerm, qs) = BilTerm(op.test_id,op.trial_id,op.test_ops,op.trial_ops,op.coeff,assign_quadstrat(op.kernel,qs))

# assign_quadstrat(op::BilForm, qs, types) = BilForm(op.test_space,op.trial_space,assign_quadstrat(op.terms,qs,types))
# assign_quadstrat(op::BilTerm, qs, types) = BilTerm(op.test_id,op.trial_id,op.test_ops,op.trial_ops,op.coeff,assign_quadstrat(op.kernel,qs,types))
# scalartype(a::AssignQuadStrat) = scalartype(a.operator_types)

function assemble!(op::AssignQuadStrat,test_functions::Space, trial_functions::Space,
    store, threading::Type{Threading{:multi}};quadstrat=nothing)
    assemble!(op.op,test_functions,trial_functions,store,threading; quadstrat=op.qs)
end

function assemble(op::AssignQuadStrat, test_functions, trial_functions;
    storage_policy = Val{:bandedstorage},
    threading = Threading{:multi},
    quadstrat=nothing)
    assemble(op, test_functions, trial_functions; storage_policy=storage_policy, threading=threading, quadstrat=op.qs)
end