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
    return exp(log(fitness[1])*p + log(fitness[2])*p)
end

function geometrical_average(fitness;p=0.5)
    return (fitness[1]*fitness[2])^(p)
end

gen_fitness = (dry = 0.785, wet = 0.785)
dry_fitness = (dry = 1.0, wet = 0.58)
wet_fitness = (dry = 0.6, wet = 1.0)

