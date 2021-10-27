using Symbolics, SymbolicUtils, LinearAlgebra, BlockDiagonals

@syms h[..] hc[..] v[..] vc[..]

"""
    DensityOperator(state)

Converts a pure state given by a polynomial in creation operators into a density operator given by a polynomial in creation and annihilation operators.
"""
function DensityOperator(state)
    conjState = state

    for n = 1:10
        conjState = substitute(conjState, Dict(h[n] => hc[n], v[n] => vc[n]))
    end

    rho = expand(state*conjState)
    return rho
end

"""
    ApplyU(state, U, modes)

Apply unitary U to state on a specified set of modes.
"""
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

"""
    ApplyU(state, U)

Apply unitary U to state on the top set of modes.
"""
function ApplyU(state, U)
    dim = Int(size(U)[1]/2)
    return ApplyU(state, U, 1:dim)
end

"""
    PhaseShifter(phase)

Returns unitary matrix for a phase-shifter.
"""
function PhaseShifter(phase)
    M = [exp(phase*im) 0; 0 exp(phase*im)]
    return M
end

"""
    PhaseShifter(phaseH,phaseV)

Returns unitary matrix for a phase-shifter with polarization-dependent phases phaseH and phaseV.
"""
function PhaseShifter(phaseH, phaseV)
    M = [exp(phaseH*im) 0; 0 exp(phaseV*im)]
    return M
end

"""
    BeamSplitter(eta)

Returns unitary matrix for a beamsplitter with reflectivity eta.
"""
function BeamSplitter(eta)
    M = [sqrt(eta) sqrt(1-eta); sqrt(1-eta) -sqrt(eta)]
    return BlockDiagonal([M,M])
end

"""
    PBS()

Returns unitary matrix for a polarising beam splitter with horizontal and vertical reflectivities etaH=1 and etaV=0.
"""
function PBS()
    H = [1 0; 0 1]
    V = [0 1; 1 0]
    return BlockDiagonal([H,V])
end

"""
    PBS(etaH,etaV)

Returns unitary matrix for a polarising beam splitter with horizontal and vertical reflectivities etaH and etaV.
"""
function PBS(etaH, etaV)
    H = [sqrt(etaH) sqrt(1-etaH); sqrt(1-etaH) -sqrt(etaH)]
    V = [sqrt(etaV) sqrt(1-etaV); sqrt(1-etaV) -sqrt(etaV)]
    return BlockDiagonal([H,V])
end

"""
    Rotate(theta)

Returns the unitary matrix for a polarization rotation given by angle theta.
"""
function Rotate(theta)
    M = [cos(theta) im*sin(theta); im*sin(theta) cos(theta)]
    return M
end

"""
    Hadamard()

Returns unitary matrix for a Hadamard polarization rotation.
"""
function Hadamard()
    return Rotate(pi/4)
end

"""
    Flip()

Returns unitary matrix for a polarization flip.
"""
function Flip()
    M = [0 1; 1 0]
    return M
end

"""
    Project(state, mode, power, HorV)

Project a given spatial and polarization mode from a state onto a photon-number given by power.
"""
function Project(state, mode, power, HorV='h')
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

# state = DensityOperator(h[1])
# state = ApplyU(state, BeamSplitter(0.5))
# println(state)
# state = ApplyU(state, BeamSplitter(0.5))
# state = ApplyU(state, PhaseShifter(1.0))
# println(state)
