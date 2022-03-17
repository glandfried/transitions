# @article{frank1990-variableEnvironment,
#   title={Evolution in a variable environment},
#   author={Frank, Steven A and Slatkin, Montgomery},
#   journal={The American Naturalist},
#   volume={136},
#   number={2},
#   pages={244--260},
#   year={1990},
#   publisher={University of Chicago Press},
#   url={https://www.researchgate.net/profile/Steven-Frank-4/publication/242389945_Evolution_in_a_Variable_Environment/links/597a82370f7e9b0469bf20c4/Evolution-in-a-Variable-Environment.pdf}
# }

using Test
using Random
using Statistics
using Distributions
using LinearAlgebra

pE = 0.71 # probabilidad del ambiente

A0 = (11, 10)
A1 = (15, 6)
A2 = (19, 2)

function environment(p=pE,n=1)
    rand(Binomial(n,p))
end
function fitness(r, allel)
    return r==0 ? allel[1] : allel[2]
end

function expected_value(states, p)
    return p*states[1] + (1-p)*states[2]
end
function deviation(allel, p)
    mu = expected_value(allel, p)
    return (allel[1] - mu, allel[2] - mu)
end
function variance(state, p)
    return expected_value(state.^2, p) - expected_value(state, p).^2
end

# Dymanic
function sum_deviation(A, n, p = pE)
    d = deviation(A, p)
    success = environment(p, n)
    return success*d[1] + (n-success)*d[2]
end
function next_size(A, n, p = pE)
    return n*expected_value(A, p) + sum_deviation(A, n, p)
end
function r(A, n, p = pE )
    next_size(A,n,p)/n
end


n0 = 1 # n haploid allel A0
n1 = 1 # n haploid allel A1
n2 = 1 # n haploid allel A2

n0 = next_size(A0, n0)
n1 = next_size(A1, n1)
n2 = next_size(A2, n2)
N = n0 + n1 + n2
q0, q1, q2 = n0/N, n1/N, n2/N


mu0 = expected_value(A0, pE)
mu1 = expected_value(A1, pE)
mu2 = expected_value(A2, pE)
dv0 = variance(deviation(A0, pE), pE)
dv1 = variance(deviation(A1, pE), pE)
dv2 = variance(deviation(A2, pE), pE)


expected_value(A0, pE)/expected_value(A0, pE)



