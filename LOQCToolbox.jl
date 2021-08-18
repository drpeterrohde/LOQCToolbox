using Symbolics, SymbolicUtils, LinearAlgebra, BlockDiagonals

@syms h[..] hc[..] v[..] vc[..]

function DensityOperator(state)
    conjState = state

    for n = 1:10
        conjState = substitute(conjState, Dict(h[n] => hc[n], v[n] => vc[n]))
    end

    rho = expand(state*conjState)
    return rho
end

function ApplyU(state, U, modes)
    dim = Int(size(U)[1]/2)
    rules = Dict()

    for i = 1:dim
        h_rhs = 0
        hc_rhs = 0
        v_rhs = 0
        vc_rhs = 0

        for j = 1:dim
            h_rhs += U[i,j]*h[modes[j]] + U[i,j+dim]*v[modes[j]]
            hc_rhs += conj(U[i,j])*hc[modes[j]] + conj(U[i,j+dim])*vc[modes[j]]
            v_rhs += U[i+dim,j]*h[modes[j]] + U[i+dim,j+dim]*v[modes[j]]
            vc_rhs += conj(U[i+dim,j])*hc[modes[j]] + conj(U[i+dim,j+dim])*vc[modes[j]]
        end

        push!(rules, h[modes[i]] => h_rhs)
        push!(rules, hc[modes[i]] => hc_rhs)
        push!(rules, v[modes[i]] => v_rhs)
        push!(rules, vc[modes[i]] => vc_rhs)
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

function BeamSplitter(eta)
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

function ForcedMeasure(state, mode, power, HorV='h')
    if power == 0
        if HorV == 'h'
            meas = expand(substitute(state, Dict(h[mode] => 0)))
        elseif HorV == 'v'
            meas = expand(substitute(state, Dict(v[mode] => 0)))
        end
    elseif power == 1
        sub = state
        for i = 2:10
            if HorV == 'h'
                sub = substitute(sub, Dict(h[mode]^i => 0))
            elseif HorV == 'v'
                sub = substitute(sub, Dict(v[mode]^i => 0))
            end
        end
        if HorV == 'h'
            diff = sub - substitute(sub, Dict(symbol => 0))
            red = substitute(diff, Dict(symbol => 1))
        elseif HorV == 'v'
        end
        meas = expand(red * sqrt(factorial(power)))
    else
        error("Not implemented")
    end

    return meas
end

### EXAMPLE

state = DensityOperator(h[1])
state = ApplyU(state, BeamSplitter(0.5))
println(state)
# state = ApplyU(state, BeamSplitter(0.5))
# state = ApplyU(state, PhaseShifter(1.0))
# println(state)
