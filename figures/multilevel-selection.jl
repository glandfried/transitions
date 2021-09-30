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

@doc "e: estrategia, c: cantidad de cooperadores, r: exitos, N: poblacion total"
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
function posterior_C(e, a, c, N, prior)
    r = coin(a,c)
    return coop_fitness.(e,c,r,N).*prior
end
function posterior_D(e, a, priorD, priorC)
    r = coin(a)
    return (priorC .+ priorD).*coop_fitness.(e,1,r,1)
end
function priorsCD(e,c,N)
    if c == 0
        priorsC = [(1.0.-e).+e].*0.0;
        priorsD = [(1.0.-e).+e]
    elseif (c > 0) & (c != N)
        priorsC = [(1.0.-e).+e] ./2
        priorsD = [(1.0.-e).+e] ./2
    elseif c == N
        priorsC = [(1.0.-e).+e]
        priorsD = [(1.0.-e).+e].*0.0
    end
    return (priorsC, priorsD)
end
function posterior_evidence_level_1(e,c,N,T=1000)
    # e: estrategias
    # c: cantidad de cooperadores
    # N: población total
    # T: tiempo total
    priorsC, priorsD = priorsCD(e,c,N)
    evidence = []
    for i in 1:T
        posteriorC = posterior_C(e,A,c,N, priorsC[end])
        posteriorD = posterior_D(e,A,priorsD[end], priorsC[end])
        push!(evidence,sum(posteriorC*de) + sum(posteriorD*de))
        push!(priorsC, posteriorC./evidence[end] )
        push!(priorsD, posteriorD./evidence[end] )
    end
    return (priorsC,priorsD,evidence)
end
function posterior_level_2(e,NN = 10,T=100)
    # NN: hasta que tamaño de porblación
    posteriorsM = []
    for N in 2:NN
        push!(posteriorsM, []) 
        for c in 0:N
            priorsC, priorsD, evidence = posterior_evidence_level_1(e,c,N,T)
            push!(posteriorsM[end] , prod(evidence) )
        end
    end
    return posteriorsM
end

postC, postD, predictions = posterior_evidence_level_1(e,4,5)
plot([e;e.+1.0],[postC[end];postD[end]])
plot([e;e.+1.0],[postC[10];postD[10]])
plot( predictions )

postL2 = posterior_level_2(e,10,250)
pData = sum([sum(ps) for ps in postL2])
pCoop = sum([ps[end] for ps in postL2])
pCoop/pData

function posterior_level_2_slide_coop(e, NN = 16,T=1000)
    # NN: hasta que tamaño de porblación
    posteriorsM = []
    for N in 1:NN
        priorsC, priorsD, evidence = posterior_evidence_level_1(e,N,N,T)
        push!(posteriorsM , prod(evidence) )
    end
    return posteriorsM./sum(posteriorsM)
end

postL2slide = log.(posterior_level_2_slide_coop(e,24))
plot(postL2slide ,label=false)

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



