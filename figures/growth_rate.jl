using Random
using Plots
using Statistics 

function w()
    r = rand(0:1)
    return (1.5^r)*(0.6^(1-r))
end

res = []
for _ in 1:1000
    T = 1000
    s1 = 1.0
    s2 = 1.0
    for i in 1:T
        s1 = w()*s1; s2 = w()*s2
        s1 = (s1+s2)/2
        s2 = s1
    end
    push!(res,s1)
end

exp(log(median(res))/1000)

C = []
D = []
for _ in 1:1000
    T = 500
    s1 = 1.0
    s2 = 1.0
    for i in 1:T
        s1 = w()*s1; s2 = w()*s2
        s1 = (s1)/2
        s2 += s1
    end
    push!(D,s2)
    push!(C,s1)
end

exp(log(median(D))/500)
exp(log(median(C))/500)


D = []
for _ in 1:1000
    T = 500
    s1 = 1.0
    s2 = 1.0
    for i in 1:T
        s1 = w()*s1; s2 = w()*s2
    end
    push!(D,s2)
end

exp(log(median(D))/500)
