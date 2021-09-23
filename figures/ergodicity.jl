using Plots
using Random
using Statistics
using Distributions
using LinearAlgebra

function r(p=0.5)
    rand(Binomial(1,p))
end

function f(r)
    return r==0 ? 1.5 : 0.6
end

ensamble_average = 1.5*0.5+0.6*0.5 
time_average = exp(log(1.5)*0.5+log(0.6)*0.5)

N = 33
T = 2000

e = zeros((T,N))
for i in 1:N
    e[:,i] .= cumsum(log.([ f(r()) for i in 1:T]))
end

fig = plot(e, label=false, thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="log(Recursos) ")
plot!([1,2000],[0,T*log(ensamble_average )], color="black", label=false)
#plot!([1,2000],[0,T*log(time_average )], color="black", label=false)
savefig(fig, "pdf/ergodicity_individual_trayectories.pdf")


N = 1
T = 1000000

e = zeros((T,N))
for i in 1:N
    e[:,i] .= cumsum(log.([ f(r()) for i in 1:T]))
end

fig = plot([1,T],[0,T*log(ensamble_average )], color="black", label=false, thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="log(Recursos) ")
plot!([1,T],[0,T*log(time_average )], color=1, label=false)
#plot!(e, label=false, color=1)
savefig(fig, "pdf/ergodicity_individual_trayectories_longrun.pdf")


Random.seed!(2)
N = 10000
T = 10

e = zeros((T+1,N))
for i in 1:N
    e[:,i] .= [[1]; cumprod([ f(r()) for i in 1:T])]
end

fig = plot(log.(e[:,1]), thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="log(Recursos) ", label="10^0", legend=(0.15,0.9))
plot!(1:T+1,log.([mean(e[t,1:(10^1)]) for t in 1:T+1]), linewidth=1.2, label="10^1")
plot!(1:T+1,log.([mean(e[t,1:(10^2)]) for t in 1:T+1]),linewidth=1.4, label="10^2")
plot!(1:T+1,log.([mean(e[t,1:(10^3)]) for t in 1:T+1]), linewidth=1.6, label="10^3")
plot!(1:T+1,log.([mean(e[t,1:(10^4)]) for t in 1:T+1]), label="10^4",linewidth=1.8, color="black")
savefig(fig, "pdf/ergodicity_expectedValue.pdf")

