using Plots
using Random
using Statistics
using Distributions
using LinearAlgebra

function coin(p=0.71, n=1)
    rand(Binomial(n,p))
end


function game(n=100,d=1,t=1000, seed=1; costo = 0.0, reproduccion = 0.5, muerte = 0.4)#evolutivo=true
    Random.seed!(seed)
    res = zeros((n,t+1))
    res[:,1] .= 1.0
    i=2
    while i <= t+1
        cpr = (sum(res[1:(n-d),i-1].*(1-costo))/n)        
        for a in 1:n#a=1
            r = rand([0,1])
            rate = r==0 ? (1+reproduccion) : (1-muerte)
            res[a,i] = a<=(n-d) ? cpr*rate : (cpr+res[a,i-1])*rate
        end
        i = i +1 
    end
    return res
end


A = 0.71 # ambiente
N = 10 # 
ns = [i for i in 1:N] 
de = 0.0001 # delta e (estrategia)
e = [i for i in 0.0:de:1.0] # estrategias (a veces ambiente)
"e: estrategia, c: cantidad de cooperadores, r: exitos, N: poblacion total"
function coop_fitness(e,c,r,N)
    # e: estrategia
    # c: cantidad de cooperadores
    # r: exitos
    # N: poblacion total
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
function posterior_C(e, a, c, N, prior)#a=A
    r = coin(a,c)
    return coop_fitness.(e,c,r,N).*prior
end
function posterior_D(e, a, priorD, priorC)
    r = coin(a)
    return (priorC .+ priorD).*coop_fitness.(e,1,r,1)
end
function priorG(N)
    pdf(Binomial(N,0.5))
end
function priorsCD(e,c,N)
    if c == 0
        priorsC = [(1.0.-e).+e].*0.0;
        priorsD = [(1.0.-e).+e]
    elseif (c > 0) & (c != N)
        priorsC = [(1.0.-e).+e] .* (c/N)
        priorsD = [(1.0.-e).+e] .* (1-c/N)
    elseif c == N
        priorsC = [(1.0.-e).+e]
        priorsD = [(1.0.-e).+e].*0.0
    end
    return priorsC, priorsD
end
function posterior_evidence_level_1(e,c,N,T=100, a = A)#c=9;N=9;T=100
    # e: estrategias
    # c: cantidad de cooperadores
    # N: población total
    # T: tiempo total
    pC, pD = priorsCD(e,c,N)
    priorsC = pC[1]
    priorsD = pD[1]
    joint_log_evidence = 0.0
    for _ in 1:T
        posteriorC = (c==0) ? pC[1] : posterior_C(e,a,c,N, priorsC)
        posteriorD = (c==N) ? pD[1] : posterior_D(e,a,priorsD, priorsC)
        evidence = sum(posteriorC*de) + sum(posteriorD*de)
        priorsC = posteriorC./evidence
        priorsD = posteriorD./evidence
        joint_log_evidence += log(evidence)
    end
    return priorsC,priorsD,joint_log_evidence
end
function posterior_level_2(e,NN = 10,T=100)
    # NN: hasta que tamaño de porblación
    posteriorsM = []
    for N in 2:NN
        push!(posteriorsM, []) 
        for c in 0:N
            priorsC, priorsD, log_evidence = posterior_evidence_level_1(e,c,N,T)
            push!(posteriorsM[end] , log_evidence )
        end
        posteriorsM[end].+log.(priorG(N))
    end
    return posteriorsM
end

function omega_desertor(f_c, f_d, t)
    return f_d^(t) + sum([ (f_c^i)*prod([f_d for j in (i+1):t]) for i in 1:(t)])
end



postC, postD, joint_log_evidence = posterior_evidence_level_1(e,9,9,10000)
fig = plot(e,postC, thickness_scaling = 2, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy", color=3, linewidth=2)
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion", color=1, linewidth=2)
savefig(fig, "pdf/multilevel-selection-1.pdf")
savefig(fig, "png/multilevel-selection-1.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-1.pdf pdf/multilevel-selection-1.pdf`) 

Random.seed!(2)
postC, postD, joint_log_evidence = posterior_evidence_level_1(e,8,9,10000)
fig = plot(e,postC, thickness_scaling = 2, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy", color=3, linewidth=2)
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion", color=1, linewidth=2)
savefig(fig, "pdf/multilevel-selection-2.pdf")
savefig(fig, "png/multilevel-selection-2.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-2.pdf pdf/multilevel-selection-2.pdf`) 

postC, postD, joint_log_evidence = posterior_evidence_level_1(e,7,9,10000)
fig = plot(e,postC, thickness_scaling = 2, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy", color=3, linewidth=2)
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion", color=1, linewidth=2)
savefig(fig, "pdf/multilevel-selection-3.pdf")
savefig(fig, "png/multilevel-selection-3.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-3.pdf pdf/multilevel-selection-3.pdf`) 
# P(poblaciones que contengan desertores) 
postL2 = posterior_level_2(e,10,100)
pData = sum([sum(exp.(le)) for le in postL2])
pCoop = sum([exp(le[end]) for le in postL2])
pCoop/pData

function posterior_level_2_slide_coop(e, NN = 16,T=1000)
    # NN: hasta que tamaño de población
    log_posteriorsM = []
    for N in 1:NN
        _, _, log_evidence = posterior_evidence_level_1(e,N,N,T)
        push!(log_posteriorsM , log_evidence )
    end
    return log_posteriorsM
end

postL2slide_maximum_temporal = [ maximum(coop_temporal_average.(i, e, A, i)) for i in 1:16]
postL2slide = posterior_level_2_slide_coop(e,16,10000)
fig = plot(1:16, exp.(postL2slide./10000 ), label="gmean(Evidencia)", thickness_scaling = 1.5, grid=false, ylab="Tasa de crecimiento", xlab="Tamaño de la población", legend=:bottomright, foreground_color_legend = nothing)
plot!(1:16, postL2slide_maximum_temporal, label="max(Analítica)")
savefig(fig, "pdf/multilevel-selection-4.pdf")
savefig(fig, "png/multilevel-selection-4.png")

######################


f_C = coop_temporal_average(100, 1.5/2.1, 0.5, 100)*2.1
f_c = coop_temporal_average(99, 1.5/2.1, 0.5, 100)*2.1
f_d = coop_temporal_average(1, 1.5/2.1, 0.5, 1)*2.1
T = 1000
trayectorias = log.(transpose(game(100,1,T)))
trayectorias_C = log.(transpose(game(100,0,T)))
fig= plot(trayectorias[:,1], label=false, color=3,legend=:best,foreground_color_legend = nothing, ylab="log recursos", xlab="Tiempo", linewidth=1, grid=false,thickness_scaling = 1.5)
#plot!(trayectorias[:,1], label=false, color=1, linewidth=1.5)
plot!(trayectorias[:,end], label=false, color=1, linewidth=1.5)
plot!(log.(omega_desertor.(f_c, f_d, 0:(T))), label=false,color="black",linewidth=1.5)
plot!( log.(f_c.^[i for i in 0:T]), label=false, color="black", linewidth=1.5)
plot!( log.(f_C.^[i for i in 0:T]), label=false, color="black", linestyles=:dash)
savefig(fig, "pdf/multilevel-selection-5.pdf")
savefig(fig, "png/multilevel-selection-5.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-5.pdf pdf/multilevel-selection-5.pdf`) 


######################
# Bayesian inference (multilevel biomass proportion)

# Biomasa de cada individuo por region
b_eg0 = game(2,0,150,1).*(1/4*1/2) # enteramente cooperadora
b_eg1 = game(2,1,150,1).*(1/2*1/2)
b_eg2 = game(2,2,150,1).*(1/4*1/2)

# Biomasa por grupo (posterior nivel 2)
b_g0 = [ sum(c) for c in eachcol(b_eg0)]
b_g1 = [ sum(c) for c in eachcol(b_eg1)]
b_g2 = [ sum(c) for c in eachcol(b_eg2)]
b_g = [transpose(b_g0); transpose(b_g1); transpose(b_g2)]

# Biomasa cooperadores (posterior multilevel)
w_coop = [ sum(c) for c in eachcol(b_eg0)]
w_coop  = w_coop .+ b_eg1[1,:]

# Biomasa total
B = [ sum(c) for c in eachcol(b_g)]

# P(g|a)
p = plot(b_g0./B, label="CC", thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="P( g | a1, ..., at)", color = 3, foreground_color_legend = nothing)
plot!(b_g1./B, label="CD", color=2)
plot!(b_g2./B, label="DD", color=1)

savefig(p, "png/multilevel-selection-6.png") 
savefig(p, "pdf/multilevel-selection-6.pdf") 

# P(coop|a)
p = plot(w_coop./B, label="Cooperator", thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="P( coop | a1, ..., at)", color = 3, foreground_color_legend = nothing)

savefig(p, "png/multilevel-selection-multilevel-posterior.png") 
savefig(p, "pdf/multilevel-selection-multilevel-posterior.pdf") 


# Biomasa individuos engrupo mixto (posterior nivel 1)

p = plot(b_eg1[1,:]./b_g1, label="Cooperator", thickness_scaling = 1.5, grid=false, xlab="Tiempo", ylab="P( i | a1, ..., at)", color = 3, foreground_color_legend = nothing)
plot!(b_eg1[2,:]./b_g1, label="Desertor", color=2)
savefig(p, "png/multilevel-selection-level-1-posterior.png") 
savefig(p, "pdf/multilevel-selection-level-1-posterior.pdf") 



#######################

fD = []
fC = []
N= 1000
nx = 1:N
t = 100
f_d = coop_temporal_average(1, 1.5/2.1, 0.5, 1)
for n in nx#n=950
    push!(fC,coop_temporal_average(n, 1.5/2.1, 0.5, N))
    wD101 = omega_desertor(fC[end],f_d,t+1)
    wD100 = omega_desertor(fC[end],f_d,t)
    push!(fD,wD101/wD100)
end

fig=plot((nx/N), reverse(fC),label="Cooperador", color=3, legend=(0.2,0.3),foreground_color_legend = nothing, ylab="Fitness", xlab="Proporción desertores (total)", linewidth=1.5, grid=false,thickness_scaling = 1.5)
plot!(nx/N,reverse(fD), label="Desertor", color=1)
savefig(fig, "pdf/multilevel-selection-7.pdf")
savefig(fig, "png/multilevel-selection-7.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-7.pdf pdf/multilevel-selection-7.pdf`) 

fD = []
fC = []
N= 2
nx = 0:N
t = 100
f_d = coop_temporal_average(1, 1.5/2.1, 0.5, 1)
for n in nx#n=950
    push!(fC,coop_temporal_average(n, 1.5/2.1, 0.5, N))
    wD101 = omega_desertor(fC[end],f_d,t+1)
    wD100 = omega_desertor(fC[end],f_d,t)
    push!(fD,wD101/wD100)
end

fig=scatter([(nx/(N-1))[1:end-1];(nx/(N-1))[1:end-1]], [reverse(fC)[1:end-1];reverse(fD)[2:end]], color=[3,3,1,1],foreground_color_legend = nothing, ylab="Fitness", xlab="Proporción desertores (otros)", linewidth=1.5,label=false, grid=false,thickness_scaling = 1.5 )
savefig(fig, "pdf/multilevel-selection-8.pdf")
savefig(fig, "png/multilevel-selection-8.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-8.pdf pdf/multilevel-selection-8.pdf`) 


fD = []
fC = []
N= 16
nx = 0:N
t = 100
f_d = coop_temporal_average(1, 1.5/2.1, 0.5, 1)
for n in nx#n=950
    push!(fC,coop_temporal_average(n, 1.5/2.1, 0.5, N))
    wD101 = omega_desertor(fC[end],f_d,t+1)
    wD100 = omega_desertor(fC[end],f_d,t)
    push!(fD,wD101/wD100)
end

fig=scatter((nx/(N-1))[1:end-1], reverse(fC)[1:end-1], color=3,foreground_color_legend = nothing, ylab="Fitness", xlab="Proporción desertores (otros)", linewidth=1.5,label=false, grid=false,thickness_scaling = 1.5, xlim=[0,0.21], ylim=[0.4,0.5])
scatter!((nx/(N-1))[1:end-1], reverse(fD)[2:end], color=1, label=false)
savefig(fig, "pdf/multilevel-selection-9.pdf")
savefig(fig, "png/multilevel-selection-9.png")
run(`pdfcrop --margins '0 0 0 0' pdf/multilevel-selection-9.pdf pdf/multilevel-selection-9.pdf`) 


#########

# De nuevo la tasa de crecimiento del desertor mixto

function g_DC(gDD,gCD,T)
    res = gDD^T 
    for t in 1:(T-1)
        res += gCD^t * gDD^(T-t-1)
    end
    if T>1
        res += gCD^T
    end
    return res
end

fig = plot([g_DC(0.45,0.48/2,T)/g_DC(0.45,0.48/2,T-1) for T in 1:10 ])


