using Symbolics, SymbolicUtils, LinearAlgebra, BlockDiagonals

@syms h[..] v[..]

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

function PhaseShifter(phaseH, phaseV)
    M = [exp(phaseH*im) 0; 0 exp(phaseV*im)]
    return M
end

function Beamsplitter(eta)
    M = [sqrt(eta) sqrt(1-eta); sqrt(1-eta) -sqrt(eta)]
    return BlockDiagonal([M,M])
end

function PBS()
    H = [1 0; 0 1]
    V = [0 1; 1 0]
    return BlockDiagonal([H,V])
end

function PBS(etaH, etaV)
    H = [sqrt(etaH) sqrt(1-etaH); sqrt(1-etaH) -sqrt(etaH)]
    V = [sqrt(etaV) sqrt(1-etaV); sqrt(1-etaV) -sqrt(etaV)]
    return BlockDiagonal([H,V])
end

function Rotate(theta)
    M = [cos(theta) im*sin(theta); im*sin(theta) cos(theta)]
    return M
end

function Hadamard()
    return Rotate(pi/4)
end

function Flip()
    M = [0 1; 1 0]
    return M
end

### EXAMPLE

state = h[2]*v[1]
state = ApplyU(state, PBS())
println(state)
# state = ApplyU(state, BeamSplitter(0.5))
# state = ApplyU(state, PhaseShifter(1.0))
# println(state)
