using Plots

@syms theta

state = h[1]
state = ApplyU(state, BeamSplitter(0.5), [1,2])
state = ApplyU(state, PhaseShifter(theta), [1])
state = ApplyU(state, BeamSplitter(0.5), [1,2])
measured = ForcedMeasure(state,1,0)
P(theta) = abs(substitute(measured, Dict(h[2] => 1)))^2