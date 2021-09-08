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


test = true
if test
    growth_rate = (0.71^0.71)*(0.29^0.29)
    log_growth_rate = log(0.71)*0.71 + log(0.29)*0.29
    growth_rate ≈ exp(log_growth_rate )
    growth_rate ≈ exp(sum(log.([[0.71 for _ in 1:71] ; [0.29 for _ in 1:29]]))/100)
end
fig = plot(ps, (0.5.^ps).*(0.5.^(1.0.-ps)),label="0.50", xlab="Ambiente", ylab="Tasa de supervivencia", legend=(0.15,0.9), thickness_scaling = 1.5)
plot!(ps, (0.71.^ps).*(0.29.^(1.0.-ps)),label="0.71")
plot!(ps, (0.99.^ps).*(0.01.^(1.0.-ps)),label="0.99")
 # <\Delta x>
plot!(ps, 0.71.*ps.+0.29*(1.0.-ps), label=false, line=:dot, color=2) 
plot!(ps, 0.99.*ps.+0.01*(1.0.-ps), label=false, line=:dot, color=3)
# Analisis
plot!([0.5,0.5], [(0.71^0.5)*(0.29^0.5),0.71*0.5+0.29*0.5], color="black", line=:arrow, label=false)
scatter!([0.71],[0.71^0.71*0.29^0.29], color=2,label=false)
savefig(fig, "coin4.pdf")
savefig(fig, "coin4.png")

###
# Growth rate de grupos de 2 jugadores

r_CC = (((1.5+1.5)/2)^0.25)*(((1.5+0.6)/2)^0.5)*(((0.6+0.6)/2)^0.25)
log_r_CC= log((1.5+1.5)/2)*0.25 + log((1.5+0.6)/2)*0.5 + log((0.6+0.6)/2)*0.25
r_CC ≈ exp(log_r_CC)
