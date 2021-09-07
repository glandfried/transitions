using Random
using Plots
using Distributions

p = 1.5/2.1
n = 10
N = n^7

exitos = rand(Binomial(n,p))
EXITOS = rand(Binomial(N,p))

ps = 0.0:0.0001:1.0
fig = plot(ps, pdf.(Beta(1,1),ps), thickness_scaling = 2.0, legend=false, ylim=(0.0,1.1))
savefig(fig, "coin1.pdf") 
fig = plot(ps, pdf.(Beta(1+exitos,1+n-exitos),ps), thickness_scaling = 2.0, legend=false)
savefig(fig, "coin2.pdf") 
fig = plot(ps, pdf.(Beta(1+EXITOS,1+N-EXITOS),ps), thickness_scaling = 2.0, legend=false)
savefig(fig, "coin3.pdf")

