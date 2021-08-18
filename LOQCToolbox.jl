using Symbolics, SymbolicUtils, LinearAlgebra, BlockDiagonals

modes = 5
@variables h[1:modes], v[1:modes]

function ApplyU(state, U, modes)
    dim = Int(size(U)[1]/2)
    rules = Dict()

    for i = 1:dim
        h_rhs = 0
        v_rhs = 0
        for j = 1:dim
            h_rhs += U[i,j]*h[modes[j]] + U[i,j+dim]*v[modes[j]]
            v_rhs += U[i+dim,j]*h[modes[j]] + U[i+dim,j+dim]*v[modes[j]]
        end
        push!(rules, h[modes[i]] => h_rhs)
        push!(rules, v[modes[i]] => v_rhs)
    end

    newState = expand(substitute(state, rules))
    return newState
end

function ApplyU(state, U)
    dim = Int(size(U)[1]/2)
    return ApplyU(state, U, 1:dim)
end

function PhaseShifter(phase)
    M = [exp(phase*im) 0; 0 exp(phase*im)]
    return M
end

function BeamSplitter(eta)
    M = [sqrt(eta) sqrt(1-eta); sqrt(1-eta) -sqrt(eta)]
    return BlockDiagonal([M,M])
end

###

state = h[1]*h[2]
state = ApplyU(state, BeamSplitter(0.5))
#state = ApplyU(state, PhaseShifter(1.0))
println(state)
