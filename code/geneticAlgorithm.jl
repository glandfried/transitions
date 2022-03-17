using Random
using Statistics
using Distributions
using LinearAlgebra

global REGION_SIZE = 2::Int64
global TOTAL_REGIONS = 500::Int64
global CYCLE_LENGHT = 200::Int64
global STRATEGY = 0.71::Float64
global STRATEGY_MUTATION = 0.0::Float64 # Probability of change the strategy (Uniform mutation)
global BEHAVIOR_MUTATION = 0.0::Float64 # Probability of change the behavior 
global BEHAVIOR_BIAS = 0.5::Float64 # Prior probability of behaviors
global ENVIRONMENT = 0.71::Float64
global K = 2.1::Float64 # Just a constant
global GENERATIONS = 10::Int64


function coin(;p=ENVIRONMENT, n=1)
    rand(Binomial(n,p))
end
function update_individual_resources(s; q=ENVIRONMENT)
    e = coin(p=q)
    s^e*(1-s)^(1-e)
end

function init_individuals(
    region_size=REGION_SIZE,
    total_regions=TOTAL_REGIONS,
    strategy=STRATEGY,
    strategy_mutation=STRATEGY_MUTATION,
    behavior_bias=BEHAVIOR_BIAS)
    #
    strategies = Vector{Vector{Float64}}()
    behaviors =  Vector{Vector{Int64}}()
    resources = Vector{Vector{Float64}}()
    for _ in 1:total_regions
        push!(strategies, []); push!(behaviors, []); push!(resources , [])
        for _ in 1:region_size 
            push!(strategies[end],rand(Uniform()) < strategy_mutation ? rand(Uniform()) : strategy )
            push!(behaviors[end], rand(Binomial(1,behavior_bias)))
            push!(resources[end], 1.0)
        end
    end
    return strategies, behaviors, resources 
end

function reset_resources(resources, total_regions=TOTAL_REGIONS)
    for r in 1:total_regions
        resources[r].=1.0
    end
end

function life_cycle(resources, behaviors, strategies,
    region_size=REGION_SIZE,
    total_regions=TOTAL_REGIONS,
    cycle_length=CYCLE_LENGHT,
    environment=ENVIRONMENT,
    k=K)
    #
    for t in 1:cycle_length
        for r in 1:total_regions
            for i in 1:region_size 
                resources[r][i] *= k*update_individual_resources(strategies[r][i],q=environment)
            end
            if sum(behaviors[r]) > 0
                common_good = sum(resources[r].*behaviors[r])
                resources[r] = resources[r].*(1.0.-behaviors[r]) .+ common_good/region_size
            end
        end
    end
end

function reproduction(resources, behaviors, strategies,
    region_size=REGION_SIZE,
    total_regions=TOTAL_REGIONS,
    strategy_mutation=STRATEGY_MUTATION,
    behavior_mutation=BEHAVIOR_MUTATION,
    behavior_bias=BEHAVIOR_BIAS)
    #
    probabilities = vcat((resources./sum(sum(resources)))...)
    fitness = Categorical(probabilities )
    old_behaviors = vcat(behaviors...)
    old_strategies = vcat(strategies...)
    old_behaviors[rand(fitness)]
    for r in 1:total_regions#r=1
        for i in 1:region_size#i=1
            how = rand(fitness)
            behaviors[r][i] = rand(Uniform()) < behavior_mutation ? 1-old_behaviors[how] : old_behaviors[how]
            strategies[r][i] = rand(Uniform())< strategy_mutation ? rand(Uniform()) : old_strategies[how]
        end
    end
end

function genetic_algorithm(;
    region_size=REGION_SIZE,
    total_regions=TOTAL_REGIONS,
    cycle_length=CYCLE_LENGHT,
    strategy=STRATEGY,
    strategy_mutation=STRATEGY_MUTATION,
    behavior_mutation=BEHAVIOR_MUTATION,
    behavior_bias=BEHAVIOR_BIAS,
    environment=ENVIRONMENT,
    k=K,
    generations = GENERATIONS )
    
    strategies, behaviors, resources = init_individuals(region_size,
                                                        total_regions,
                                                        strategy,
                                                        strategy_mutation,
                                                        behavior_bias)
    for _ in 1:generations 
        reset_resources(resources, total_regions)
        life_cycle(resources, behaviors, strategies,
                   region_size,
                   total_regions,
                   cycle_length,
                   environment,
                   k) 
        reproduction(resources, behaviors, strategies,
                     region_size,
                     total_regions,
                     strategy_mutation,
                     behavior_mutation,
                     behavior_bias)
    end
    return behaviors, strategies
end


final_behaviors, final_strategies = genetic_algorithm(cycle_length=1000,region_size=10 )
sum(sum(final_behaviors))
