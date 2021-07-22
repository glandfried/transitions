# @article{starrfelt2012-bet,
#   title={Bet-hedging: a triple trade-off between means, variances and correlations},
#   author={Starrfelt, Jostein and Kokko, Hanna},
#   journal={Biological Reviews},
#   url={https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1469-185X.2012.00225.x?download=true}
# }

using Plots

function div_fitness(p_dry=0.44;  dry_at_wet = 0.58, wet_at_dry = 0.6 )
    return (dry_fitness = p_dry+wet_at_dry*(1-p_dry), wet_fitness = dry_at_wet*p_dry+(1-p_dry))
end

function temporal_average(fitness;p=0.5)
    return log(fitness[1])*p + log(fitness[2])*p
end

function geometrical_average(fitness;p=0.5)
    return (fitness[1]*fitness[2])^(1/2)
end

gen_fitness = (dry = 0.785, wet = 0.785)
dry_fitness = (dry = 1.0, wet = 0.58)
wet_fitness = (dry = 0.6, wet = 1.0)

plot(1.0.+temporal_average.(div_fitness.(0.01:0.01:0.99)))
plot!([0,100],[1.0+temporal_average(gen_fitness),1.0+temporal_average(gen_fitness)])
plot!(geometrical_average.(div_fitness.(0.01:0.01:0.99)))
plot!([0,100],[geometrical_average(gen_fitness),geometrical_average(gen_fitness)])



# Se parecen. Entonces las vamos a medir en la misma escala.
plot(temporal_average.(div_fitness.(0.01:0.01:0.99)) .- temporal_average(gen_fitness))
plot!((geometrical_average.(div_fitness.(0.01:0.01:0.99)) .- geometrical_average(gen_fitness)).*1.28 )

# Conclusión:
#    Parecen iguales
# Observaci'on:
#    (1.5*0.6)^(1/2) = 0.9486832980
#    1+(log(1.5)*0.5+log(0.6)*0.5) = 0.9473197421

# ToDo:
#    Verificar el par'ametro de proporcionalidad

