# @article{rice2008-stochasticPriceEquation,
#   title={A stochastic version of the Price equation reveals the interplay of deterministic and stochastic processes in evolution},
#   author={Rice, Sean H},
#   journal={BMC evolutionary biology},
#   volume={8},
#   number={1},
#   pages={1--16},
#   year={2008},
#   publisher={Springer},
#   url={https://link.springer.com/content/pdf/10.1186/1471-2148-8-262.pdf}
# }

using Test
using Random
using Statistics
using Distributions
using LinearAlgebra

pE = 0.71 # probabilidad del ambiente

A1 = (15, 6)
A2 = (19, 2)

function environment(p=pE,n=1)
    rand(Binomial(n,p))
end
function fitness(r, allel)
    return r==0 ? allel[1] : allel[2]
end

function w(allel, p = pE)
    return rand(Binomial(1,p))==0 ? allel[1] : allel[2]
end

# \phi \in \{0, 1\}
