using Plots
function rate(l,a,b,p=0.5)
    return p*log(1+a*l) + p*log(1+b*l)
end

bs = [i for i in 1.05:0.01:2.09] .-1
ls = [i for i in 0.0:0.01:1.0]
# 1.5, 0.6 es convexa
rs = rate.(ls,(1.5)-1,(0.6)-1); l = ls[argmax(rs)]
plot(rs)
# 1.05, 1.05 recto  
rs = rate.(ls,0.05,0.05)
plot(rs); l = ls[argmax(rs)]
# 1.10, 1.00 recto
rs = rate.(ls,0.1,0.0)
plot(rs); l = ls[argmax(rs)]
# 1.30, 0.8 recto
rs = rate.(ls,1.3-1,0.8-1)
plot(rs)

# Las estrategias probabilísticas tiene un laverage 0
# basados en la predicción a priori
# Pero no en el posterior ..
pd_l = rate.(ls,0.,0.)
for e in 0.01:0.01:0.99
    pd_l += rate.(ls,(e)-1,(1-e)-1,0.71).*(1/99)
end
