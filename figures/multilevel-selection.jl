using Plots
using Random
using Statistics
using Distributions
using LinearAlgebra

function coin(p=0.71, n=1)
    rand(Binomial(n,p))
end

A = 0.71 # ambiente
N = 10 # 
ns = [i for i in 1:N] 
de = 0.0001 # delta e (estrategia)
e = [i for i in 0.0:de:1.0] # estrategias (a veces ambiente)

function coop_fitness(e,c,r,N)
    return ((c-r)/N)*(1-e)+(r/N)*e
end
function coop_temporal_average(n, e, ambiente, N)
    res = 1.0
    p_fitness = pdf(Binomial(n,ambiente))
    for r in 0:n#r=0
        res *= coop_fitness(e,n,r,N)^p_fitness[r+1]
    end
    return res
end


# Posteriors cooperadores

coin(0.71,10)
c=10
r=7
N=10
function posterior_coop(e, a, c, N, prior)
    r = coin(a,c)
    return coop_fitness.(e,c,r,N).*prior
end

priors = [(1.0.-e).+e]
for i in 1:20
    posterior = posterior_coop(e,A,10,10, priors[end])
    push!(priors, posterior./sum(posterior*de) )
end
plot(e,priors[1])
plot!(e,priors[3])
plot!(e,priors[5])
plot!(e,priors[7])
plot!(e,priors[9])


# Resultado analítico, el mismo que cpr 
function omega_desertor(f_c, f_d, t)
    return sum([ f_c^i*f_d^(t-i) for i in 1:t])
end

f_c = 1.05
f_d = 1.5^0.5*0.6^0.5
T = 1000
diff = []
omega_d = []
omega_c = []
for t in 1:T 
    push!(omega_d,log(omega_desertor(f_c, f_d, t)))
    push!(omega_c, log(f_c^t))
    push!(diff,omega_d[end] - omega_c[end])
end
plot(diff)
plot!(log.(analyticTimeAverageDefectorBiomass.(0:t)) .- log.(analyticTimeAverageCoopratorBiomass.(0:t)))

plot(omega_c)
plot!(omega_d)



