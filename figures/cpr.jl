using Random
using Plots

delta_expected = 1.5*0.5+0.6*0.5 - 1 # <\Delta x>
time_average = log(1.5)*0.5+log(0.6)*0.5 # <\Delta ln x>
ensamble_average = 1 + delta_expected


ri = sqrt(1.5*0.6)
d = 0.01 
t = 100
re = 1 + delta_expected

function analyticTimeAverageCoopratorBiomass(T=t; desertor_proportion = d, ensamble_growth = re, initial_biomass = 1.0)
    return t == 0 ? initial_biomass : ((1-desertor_proportion)*ensamble_growth)^T
end

plot([0:t],log.(analyticTimeAverageCoopratorBiomass.(0:t)),legend = false)
plot!([0,t],log.([1.0,(1+delta_expected)^(t) ]), color="black")
plot!([0,t],log.([1.0,(1+time_average)^(t)]), color="black")

function analyticTimeAverageDefectorBiomass(T=t; desertor_proportion = d, individual_growth=ri, ensamble_growth = re, initial_biomass = 1.0 )
    WD0 = initial_biomass
    WCd = (1-desertor_proportion)*ensamble_growth
    ri = individual_growth
    return WD0*(ri^T) + sum([ (WCd^t)*(ri^(T-t)) for t in 1:(T-1)])
end

t=1000
p = plot([0:t],log.(analyticTimeAverageCoopratorBiomass.(0:t)),legend = false)
plot!([0:t],log.(analyticTimeAverageDefectorBiomass.(0:t)),legend = false)
plot!([0,t],log.([1.0,(1+delta_expected)^(t) ]), color="black",legend = false)
plot!([0,t],log.([1.0,(1+time_average)^(t)]), color="black",legend = false)

savefig(p, "cpr_analyticTimeAverage_d0.01.pdf") 

t=100
p = plot(analyticTimeAverageDefectorBiomass.(1:t+1)./analyticTimeAverageDefectorBiomass.(0:t))
plot!([1,t],[(1-d)*ensamble_average,(1-d)*ensamble_average])

savefig(p, "cpr_analyticTimeAverage_convergence_d0.01.pdf") 

WD = []
WC = []
dx = 1:10000
t = 1000
for i in dx#i=5000
    di = i/(dx[end]+1)
    WDd = analyticTimeAverageDefectorBiomass(t+1,desertor_proportion=di)/analyticTimeAverageDefectorBiomass(t,desertor_proportion=di)
    push!(WD,WDd )
    WCd = analyticTimeAverageCoopratorBiomass(2,desertor_proportion=di)/analyticTimeAverageCoopratorBiomass(1,desertor_proportion=di)
    push!(WC,WCd )
end

p = plot(dx./10000, WD, ylab="fitness", xlab="Proportion of defectors", legend=false)
plot!(dx./10000,WC)

savefig(p, "cpr_analyticTimeAverage_frequency-based-fitness.pdf") 

p = plot((dx./10000)[1:end-1],WD[2:end].-WD[1:end-1], legend=false, ylab="Change in fitness", xlab="Proportion of defectors")
plot!((dx./10000)[1:end-1],WC[2:end].-WC[1:end-1])

savefig(p, "cpr_analyticTimeAverage_frequency-based-fitness_change-fitness-transition.pdf") 

function timeAverageDeDesertorEnPoblacionInfinitaDeCooperadores(t=100000, pgg_growth=1.05)
    res = [0.0]
    pgg = 1.0
    for i in 1:t
        cara = rand([0,1]) == 0
        if cara
            push!(res, (res[i] + pgg)*1.5 )
        else
            push!(res, (res[i] + pgg)*0.6 )
        end
        pgg = pgg*pgg_growth
    end
    return res
end

function distribucionDelDesertorInvasor(;n=100, t=10000, pgg_growth=1.05)
    res = []
    for _ in 1:n
        push!(res, log(timeAverageDeDesertorEnPoblacionInfinitaDeCooperadores(t,pgg_growth)[end]))
    end
    return res
end

t=10000
wd = distribucionDelDesertorInvasor(n=1000,t=t,pgg_growth=1.0)
histogram(wd)

d_equilibrio = 1.0/(log(2.0)/log(1.05/0.94) - 1)

n=2;d=0;t=150; seed=1; costo = 0.0; reproduccion = 0.5; muerte = 0.4; evolutivo=true

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

#plot(log.(transpose(game(10,0,100000,reproduccion=0.05,muerte=0.04))), legend=false)
#plot(log.(transpose(game(2,2,1000,reproduccion=0.5,muerte=0.4))), legend=false)


###################
# Absolute 

p = plot((transpose(game(2,0,100))), legend=false, color="green")
plot!((transpose(game(2,1,100))), legend=false, color="blue")
plot!((transpose(game(2,2,100))), legend=false, color="red")

savefig(p, "cpr_absoulte_two.pdf") 

p = plot(log.(transpose(game(2,0,150))), legend=false, color="green")
plot!(log.(transpose(game(2,1,150))), legend=false, color="blue")
plot!(log.(transpose(game(2,2,150))), legend=false, color="red")

savefig(p, "cpr_absoulte_two_log.pdf") 

h=game(33,33,2000,9)
ylim = (log(minimum(h)),log((1+delta_expected)^2000))

p = plot(log.(transpose(h)),legend = false, thickness_scaling = 1.5, grid=false,  ylim = ylim)
plot!([1,2000],log.([1.0,(1+delta_expected)^2000 ]), color="black")
plot!([1,2000],log.([1.0,(1+time_average)^2000]), color="black")

savefig(p, "cpr_absoulte_individual.pdf") 

p = plot(log.(transpose(game(1,1,10,9))),legend = false, thickness_scaling = 1.5, grid=false)
plot!(log.([ sum(c)/length(c) for c in eachcol(game(10,10,10,9))]))
plot!(log.([ sum(c)/length(c) for c in eachcol(game(100,100,10,9))]))
plot!(log.([ sum(c)/length(c) for c in eachcol(game(1000,1000,10,9))]))
plot!(log.([ sum(c)/length(c) for c in eachcol(game(10000,10000,10,9))]), color="black")

savefig(p, "cpr_absolute_individual_media_estados.pdf") 

p = plot(log.(game(100,0,2000)[1,:]),legend = false, color="green", thickness_scaling = 1.5, grid=false, ylim = ylim,linewidth=2.5)
plot!([1,2000],log.([1.0,(1+delta_expected)^2000 ]), color="black")
plot!([1,2000],log.([1.0,(1+time_average)^2000]), color="black")

savefig(p, "cpr_absolute_cooperation.pdf") 

p = plot(log.(transpose(game(100,0,2000))),legend = false, color="green", thickness_scaling = 1.5, grid=false, ylim = ylim)
plot!(log.(transpose(game(100,1,2000))),legend = false, color="blue")
plot!(log.(transpose(game(100,10,2000))),legend = false, color="red")
#plot!([1,2000],log.([1.0,(1+delta_expected)^2000 ]), color="black")
#plot!([1,2000],log.([1.0,(1+time_average)^2000]), color="black")

savefig(p, "cpr_absolute_cooperation_defection.pdf") 

p = plot(log.(transpose(game(10,0,10,1,costo=0.01))), legend = false, color="green",thickness_scaling = 1.5, grid=false)
plot!(log.(transpose(game(10,1,10,1,costo=0.01))), legend = false, color="blue")
plot!(log.(transpose(game(10,2,10,1,costo=0.01))), legend = false, color="red")

savefig(p, "cpr_absolute_cooperation_defection_costo_zoom.pdf") 

p = plot(log.(transpose(game(100,0,2000,costo=0.01))),legend = false, color="green", thickness_scaling = 1.5, grid=false, ylim = ylim)
plot!(log.(transpose(game(100,1,2000,costo=0.01))),legend = false, color="blue")
plot!(log.(transpose(game(100,10,2000,costo=0.01))),legend = false, color="red")

savefig(p, "cpr_absolute_cooperation_defection_costo.pdf") 
savefig(p, "cpr_absolute_cooperation_defection_costo.png") 

######################
# Bayesian inference (multilevel biomass proportion)

# Biomasa de cada individuo por grupo
b_eg0 = game(2,0,10)
b_eg1 = game(2,1,10)
b_eg2 = game(2,2,10)

# Biomasa por grupo
b_g0 = [ sum(c) for c in eachcol(b_eg0)]
b_g1 = [ sum(c) for c in eachcol(b_eg1)]
b_g2 = [ sum(c) for c in eachcol(b_eg2)]
b_g = [transpose(b_g0); transpose(b_g1); transpose(b_g2)]

# Biomasa total
B = [ sum(c) for c in eachcol(b_g)]

# P(g|r)
p = plot(b_g0./B, label=false)
plot!(b_g1./B, label=false)
plot!(b_g2./B, label=false)

savefig(p, "cpr_bayesian_inference_multilevel_biomass_p(g|r).png") 
savefig(p, "cpr_bayesian_inference_multilevel_biomass_p(g|r).pdf") 

# P(e|g,r)
p = plot(b_eg1[1,:]./b_g1, label=false)
plot!(b_eg1[2,:]./b_g1, label=false)

savefig(p, "cpr_bayesian_inference_multilevel_biomass_p(e|g,r).png") 
savefig(p, "cpr_bayesian_inference_multilevel_biomass_p(e|g,r).pdf") 

# Vuelvo a escribir la función game()
# pero desde un perspectiva bayesiana

function w(r)
    return (1.5^(r))*0.6^(1-r)
end

function bayesian_inference_process(;N = 6, T = 250, intercalar = false)
    R = zeros(Int64,(2,T)) 
    R[1,:] = intercalar ? [ mod(i,2) for i in 1:T] : [ rand([0,1]) for _ in 1:T]
    R[2,:] = intercalar ? [ mod(i+1,2) for i in 1:T] : [ rand([0,1]) for _ in 1:T]
    e = zeros((N,T+1))
    e[:,1] .= 1/6
    g = zeros((3,T))
    g[:,1] .= 1/3
    for t in 2:(T)#t=2
        
        # P(et|r) COOP
        e[1,t] = (e[1,t-1]+e[2,t-1])/2
        e[2,t] = (e[1,t-1]+e[2,t-1])/2
        # P(et|r) MIX
        e[3,t] = e[3,t-1]/2
        e[4,t] = e[3,t-1]/2 + e[4,t-1]
        # P(et|r) DEFECT
        e[5,t] = e[5,t-1]
        e[6,t] = e[6,t-1]
        
        # P(e|r) \propto
        e[1,t] = e[1,t]*w(R[1,t-1])
        e[2,t] = e[2,t]*w(R[2,t-1])
        e[3,t] = e[3,t]*w(R[1,t-1])
        e[4,t] = e[4,t]*w(R[2,t-1])
        e[5,t] = e[5,t]*w(R[1,t-1])
        e[6,t] = e[6,t]*w(R[2,t-1])
        
        # P(e|r) = [normalization]
        e[:,t] = e[:,t]./sum(e[:,t])
        
        # p(g|e)
        g[1,t] = e[1,t] + e[2,t]
        g[2,t] = e[3,t] + e[4,t]
        g[3,t] = e[5,t] + e[6,t]
        
    end
    return e, g
end
e, g = bayesian_inference_process(T=100,intercalar=true)   
plot(transpose(e))
plot(transpose(g))

###################
# Relative Level 1 (interior del grupo)

function relative(g)
    res = zeros(size(g))
    denom = [ sum(c) for c in eachcol(g)]
    for i in 1:size(g)[1]
        res[i,:] .= g[i,:]./denom
    end
    return res
end

r = relative(game(10,0,200,costo=0.01))
p = plot(log.(transpose(r[end-5:end,:])),legend = false, color="green", thickness_scaling = 1.5, grid=false)
r = relative(game(10,1,200,costo=0.01))
plot!(log.(transpose(r[end-5:end,:])),legend = false, color="blue", thickness_scaling = 1.5, grid=false)
r = relative(game(10,2,200,costo=0.01))
plot!(log.(transpose(r[end-5:end,:])),legend = false, color="red", thickness_scaling = 1.5, grid=false)

savefig(p, "cpr_intra-relative_cooperation_defection_costo.pdf") 

###################
# Relative Level 2 (entre grupos)

function absoluto_global(g)
    return [sum(c) for c in eachcol(g)]
end

pops = zeros((3,201))
pops[1,:] .= absoluto_global(game(10,0,200,costo=0.01))
pops[2,:] .= absoluto_global(game(10,1,200,costo=0.01))
pops[3,:] .= absoluto_global(game(10,2,200,costo=0.01))

pops_relative = relative(pops)
p = plot(log.(pops_relative[1,:]),legend = false, thickness_scaling = 1.5, grid=false, color="green")
plot!(log.(pops_relative[2,:]), color="blue")
plot!(log.(pops_relative[3,:]), color="red")

savefig(p, "cpr_inter-relative_cooperation_defection_costo.pdf") 

#######################
# Reproductive multiplicatve process
#
# Si bien el proceso multiplicativo propuesto por Ole Peters puede interpretarse como secuencia de reproducci'on y superviviencia, el proceso de bienes comunes no considera la frecuencia relativa de las estrategias, sino la valor absoluto.
# En este ejemplo vamos a trabajar con una cantidad de tiradores de moneda (agentes) finita y constante (ej n=100): cada uno de los agentes comienzan una estrategia, por ejemplo 99% cooperadores y 1% desertor.
# En este ejemplo, luego de cada iteración, vamos a actualizar la estrategia de los agentes en función de la proporción que represeta cada clase de estrategias.

h, A, B = game(100,1,1000,evolutivo=true)

p = plot(A, legend = false, thickness_scaling = 1.5)

savefig(p, "cpr_reproduccion_proporcion_desertora.pdf") 

p = plot(log.(transpose(h)), legend = false, thickness_scaling = 1.5, grid=false)

savefig(p, "cpr_reproduccion_absoluto.pdf") 


#######################

function dilema(n=100,d=1,t=1000; b=50, c=1)
    res = zeros((n,t+1))
    res[:,1] .= 0.0
    for i in 2:(t+1)#i=2
        cpr = sum(b*(n-d))/n
        for a in 1:n#a=1
            res[a,i] = res[a,i-1] + cpr - (a<=(n-d) ? c : 0)
        end
    end
    return res
end

sum(dilema(10,0)[:,end])
sum(dilema(10,1)[:,end])
sum(dilema(10,2)[:,end])

p = plot(transpose(dilema(10,0,c=1)),legend=false, color="green", linewidth=2)
plot!(transpose(dilema(10,1,c=1)[1:(end-1),:]),legend=false, color="blue",linewidth=2)
plot!(dilema(10,1,c=25)[end,:],legend=false, color="blue")
plot!(transpose(dilema(10,2,c=1)[1:(end-2),:]),legend=false, color="red",linewidth=2)
plot!(transpose(dilema(10,2,c=1)[(end-1):end,:]),legend=false, color="red")

savefig(p, "cpr_prisioner_dilema_lowcost.pdf") 
savefig(p, "cpr_prisioner_dilema_lowcost.png") 

p = plot(transpose(dilema(10,0,c=25)),legend=false, color="green", linewidth=2)
plot!(transpose(dilema(10,1,c=25)[1:(end-1),:]),legend=false, color="blue",linewidth=2)
plot!(dilema(10,1,c=25)[end,:],legend=false, color="blue")
plot!(transpose(dilema(10,2,c=25)[1:(end-2),:]),legend=false, color="red",linewidth=2)
plot!(transpose(dilema(10,2,c=25)[(end-1):end,:]),legend=false, color="red")


savefig(p, "cpr_prisioner_dilema_highcost.pdf") 
savefig(p, "cpr_prisioner_dilema_highcost.png") 



