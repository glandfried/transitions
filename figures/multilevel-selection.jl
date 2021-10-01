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
    evidence = 1.0
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
    end
    return posteriorsM
end

postC, postD, joint_log_evidence = posterior_evidence_level_1(e,9,9,10000)
fig = plot(e,postC, thickness_scaling = 1.5, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy")
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion")
savefig(fig, "pdf/multilevel-selection-1.pdf")
savefig(fig, "png/multilevel-selection-1.png")

postC, postD, joint_log_evidence = posterior_evidence_level_1(e,8,9,10000)
fig = plot(e,postC, thickness_scaling = 1.5, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy")
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion")
savefig(fig, "pdf/multilevel-selection-2.pdf")
savefig(fig, "png/multilevel-selection-2.png")

postC, postD, joint_log_evidence = posterior_evidence_level_1(e,7,9,10000)
fig = plot(e,postC, thickness_scaling = 1.5, grid=false, label="Cooperation", legend=:best,foreground_color_legend = nothing, ylab="Density", xlab="Estrategy")
plot!(-1.0.*reverse(e),reverse(postD), label="Desertion")
savefig(fig, "pdf/multilevel-selection-3.pdf")
savefig(fig, "png/multilevel-selection-3.png")


# P(poblaciones que contengan desertores) 
postL2 = posterior_level_2(e,10,100)
pData = sum([sum(exp.(le)) for le in postL2])
pCoop = sum([exp(le[end]) for le in postL2])
pCoop/pData

function posterior_level_2_slide_coop(e, NN = 16,T=1000)
    # NN: hasta que tamaño de porblación
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

function omega_desertor(f_c, f_d, t)
    return sum([ (f_c^i)*(f_d^(t-i)) for i in 1:t])
end

function game(n=100,d=1,t=1000, seed=1; costo = 0.0, reproduccion = 0.5, muerte = 0.4, evolutivo=false, intercalar=false)#evolutivo=true
    Random.seed!(seed)
    res = zeros((n,t+1))
    res[:,1] .= 1.0
    desertores = []
    proporcion = []
    i=2
    while i <= t+1
        if evolutivo & (d != n) & (d != 0) 
            p = sum((res[:,i-1]./sum(res[:,i-1]))[1:(n-d)])
            d = convert(Int64,round((1-p)*n))
            push!(desertores,d)
            push!(proporcion,p)
        end
        cpr = (sum(res[1:(n-d),i-1].*(1-costo))/n)        
        for a in 1:n#a=1
            r = rand([0,1])
            if (intercalar & (mod(a,2) == mod(i,2))) || (!intercalar & r == 0)
                res[a,i] = a<=(n-d) ? cpr*(1+reproduccion) : (cpr+res[a,i-1])*(1+reproduccion)
            elseif (intercalar & (mod(a,2) != mod(i,2)) ) || (!intercalar & r == 1)
                res[a,i] = a<=(n-d) ? cpr*(1-muerte) : (cpr+res[a,i-1])*(1-muerte)
            end
        end
        i = i +1 
    end
    return evolutivo ? (res, desertores, proporcion) : res
end


f_c = coop_temporal_average(99, 0.71, 0.5, 100)*2.1
f_d = 1.5^0.5*0.6^0.5
T = 1000
fig= plot(log.(transpose(game())), label=false)
plot!(log.(omega_desertor.(f_c, f_d, 0:T)), label=false)
plot!( log.(f_c.^[i for i in 0:T]), label=false, color="black")
savefig(fig, "pdf/multilevel-selection-5.pdf")
savefig(fig, "png/multilevel-selection-5.png")


