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
plot!([1,2000],[0,T*log(time_average )], color="black", label=false)
savefig(fig, "pdf/ergodicity_individual_trayectories.pdf")

