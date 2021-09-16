using Plots
using Random
using Statistics
using Distributions
using LinearAlgebra

function coin(p=0.71)
    rand(Binomial(1,p))
end

A = 0.71
N = 10
ns = [i for i in 1:N]
e = [i for i in 0.0:0.0001:1.0]

function coop_fitness(e,n,r)
    return ((n-r)/n)*(1-e)+(r/n)*e
end

function coop_temporal_average(n, e, ambiente=A)
    log_res = 0.0
    p_fitness = pdf(Binomial(n,ambiente))
    for r in 0:n#r=0
        log_res += log(coop_fitness(e,n,r))*p_fitness[r+1]
    end
    return exp(log_res)
end

coop_temporal_average.(100,)

eN = zeros(Float64,(N,length(e)))
for i in 1:N
    eN[i,:] .= coop_temporal_average.(i,e)
end
plot(e,eN[1,:])
plot!(e,eN[2,:])
plot!(e,eN[3,:])
plot!(e,eN[4,:])
plot!(e,eN[5,:])
plot!(e,eN[6,:])
plot!(e,eN[10,:])


aN71 = zeros(Float64,(N,length(e)))
for i in ns
    aN71[i,:] .= coop_temporal_average.(i,0.71,e)
end
plot(e,aN71[1,:])
plot!(e,aN71[2,:])
plot!(e,aN71[3,:])
plot!(e,aN71[4,:])
plot!(e,aN71[5,:])
plot!(e,aN71[6,:])
plot!(e,aN71[7,:])
plot!(e,aN71[10,:])
plot!(e, 0.71.*e.+0.29*(1.0.-e), label=false, line=:dot, color=2)

aN = zeros(Float64,(length(n),length(estrategias)))
for i in 1:N
    aN[i,:] .= coop_temporal_average.(i,0.99,e)
end
plot(e,aN[1,:], legend=(0.15,0.9))
plot!(e,aN[2,:])
plot!(e,aN[3,:])
plot!(e,aN[4,:])
plot!(e,aN[5,:])
plot!(e,aN[6,:])
plot!(e,aN[7,:])
plot!(e,aN[8,:])
plot!(e,aN[9,:])
plot!(e,aN[10,:])
plot!(e, 0.99.*e.+0.01*(1.0.-e), label=false, line=:dot, color=2)


