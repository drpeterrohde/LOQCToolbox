include("LOQCToolbox.jl")


### State creation ###
# We define a single photon in mode m with 'H' polarisation as h[m] = [1, 0]^Transpose
# We define a single photon in mode m with 'V' polarisation as v[m] = [0, 1]^Transpose
# hc[m] and vc[m] refer to conjugates of h[m] and v[m] respectively = [1, 0] and [0, 1] respectively
# One can think of h[m] and v[m] as column vectors where hc[m] and vc[m] are their transpose conjugate row vectors
# Examples:

state1 = (h[1] + v[1])/sqrt(2) #superposition state of 1 qubit in 1 mode
state2 = (h[1]*h[2] + v[1]*v[2])/sqrt(2) #Bell-pair between mode 1 and 2
println("State: ", state2)


### Density Matrix creation ###
# Use function DensityOperator(state) to create density matrix out of a state
# Example

dm1 = DensityOperator(state1)
println("Density operator: ", dm1)


### Defining fundamental unitaries in optics ###

## 1. Phase shifters ##
PS1 = PhaseShifter(pi/4)
println("Phase-shifter matrix: ", PS1)

## 2. Phase Shifter - different phase shifts for different polarisations ##
PS2 = PhaseShifter(pi/4, pi/2)
println("General Phase-shifter matrix: ", PS2)

## 3. Beam-splitter ##
BS = BeamSplitter(0.5)
println("Beam-splitter matrix: ", BS)

## 4. Polarising Beam-Splitter ##
PBS1 = PBS()
println("PBS matrix: ", PBS1)

## 5. Phase Shifter - different reflectivities for different polarisations ##
PBS2 = PBS(0.4,0.5)
println("General PBS matrix: ", PBS2)

## 6. Rotation Matrix - does general single qubit rotations ##
R = Rotate(pi/6)
println("Rotation matrix: ", R)

## 7. Hadamard - special case of rotation matrix where theta=pi/4 ##
H = Hadamard()
println("Hadamard matrix: ", H)

## 8. Flip - also called X matrix, special case of rotation matrix where theta=0 ##
F = Flip()
println("Flip matrix: ", F)


### Applying unitaries on states ###

## 1. Applying single qubit unitaries ##
println("Initial state: ", state2)
state3 = ApplyU(state2, F, 1)
println("After applying flip on 1st mode: ", state3)

## 2. Applying two-qubit unitaries ##
println("Initial state: ", state2)
state4 = ApplyU(state2, PBS(0.4,0.7), [1,2])
println("After applying unbalanced PBS between 1st and 2nd mode: ", state4)


### Applying measurement on states ###

## 1. Single photon measurement ##
state5 = 3*h[1]*v[2]*h[3]*v[2] + h[2]*h[2]*h[2]*v[2] + 2*v[1]*v[2]*h[3]*v[2]
println("Initial state: ", state5)
state6 = Measure(state5, v[2])
println("State after a vertically polarised photon has been measured in mode 2: ", state6)

## 2. n-photon measurement
# Example: here n=2
state7 = Measure(state5, v[2], 2)
println("State after two vertically polarised photons has been measured in mode 2: ", state7)
