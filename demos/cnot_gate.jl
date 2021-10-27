@variables alpha beta gamma delta

eta1 = 1.0/3
eta2 = (3+sqrt(6))/6.0

input = (alpha*h[1] + beta*h[3]) * (gamma*h[2] + delta*h[4]) * h[5] * h[6]

s1 = ApplyU(input, BeamSplitter(eta1), [3,4])
s2 = ApplyU(s1, BeamSplitter(eta1), [5,6])
# s3 = ApplyU(s2, BeamSplitter(eta1), [3,4])
# s4 = ApplyU(s3, BeamSplitter(eta2), [5,6])