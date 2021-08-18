using Symbolics, LinearAlgebra, BlockDiagonals

modes = 5
@variables h[1:modes], v[1:modes]

function ApplyU(state, U::Matrix{ComplexF64}, modes)
    dim = Int(size(U)[1]/2)
    rules = Dict()

    for i = 1:dim
        h_rhs = 0.0+0.0im
        v_rhs = 0.0+0.0im
        for j = 1:dim
            h_rhs += U[i,j]*h[modes[j]] + U[i,j+dim]*v[modes[j]]
            v_rhs += U[i+dim,j]*h[modes[j]] + U[i+dim,j+dim]*v[modes[j]]
        end
        push!(rules, h[modes[i]] => h_rhs)
        push!(rules, v[modes[i]] => v_rhs)
    end

    newState = substitute(state, rules)
    return newState
end

function ApplyU(state, U::Matrix{ComplexF64})
    dim = Int(size(U)[1]/2)
    return ApplyU(state, U, 1:dim)
end

function PhaseShifter(phase::Float64)::Matrix{ComplexF64}
    M = [exp(phase*1.0im) 0.0+0.0im; 0.0+0.0im exp(phase*1.0im)]
    return M
end

function BeamSplitter(eta::Float64)::Matrix{ComplexF64}
    M = [sqrt(eta)+0.0im sqrt(1-eta)+0.0im; sqrt(1-eta)+0.0im -sqrt(eta)+0.0im]
    return BlockDiagonal([M,M])
end

###

state = h[1]*h[2]
state2 = ApplyU(state, BeamSplitter(0.5))
state3 = ApplyU(state2, PhaseShifter(1.0))
println(state3)
