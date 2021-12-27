using Random
using Plots
using Distributions

global de = 0.001
e = (de/2):de:1-(de/2)
global E = e

function f(e;a=[1])
    e^sum(a) * (1-e)^(length(a)-sum(a))
end

function pe_a(e;a)
    f.(e,a=a)./(sum(f.(E,a=a))*de )
end

pa = Binomial(1,0.71)

as = rand(pa)
plot(pe_a.(e,a=as), legend=false)
as = [as;rand(pa)] 
plot!(pe_a.(e,a=as), legend=false)
as = [as;rand(pa,8)] 
plot!(pe_a.(e,a=as), legend=false)
as = [as;rand(pa,1000)] 
plot!(pe_a.(e,a=as), legend=false)
plot!(pdf.(Beta(sum(as),length(as)-sum(as)),e)) 

# Conclusión:
#    Cuando el espacio de estrategias tiene como fitness la función
#    e^a * (1-e)^(n-a), el posterior es Beta.

############################
# Version vieja

p = 1.5/2.1
n = 10
N = 10^5

exitos = rand(Binomial(n,p))
EXITOS = rand(Binomial(N,p))


ps = 0.0:0.0001:1.0
fig = plot(ps, pdf.(Beta(1,2),ps), thickness_scaling = 2.0, legend=false)
savefig(fig, "coin1.pdf") 
fig = plot(ps, pdf.(Beta(1+exitos,1+n-exitos),ps), thickness_scaling = 2.0, legend=false)
savefig(fig, "coin2.pdf") 
fig = plot(ps, pdf.(Beta(1+EXITOS,1+N-EXITOS),ps), thickness_scaling = 2.0, legend=false)
savefig(fig, "coin3.pdf")
