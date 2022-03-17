# Factor graph     Graphical model
# [ ]       [ ] P(a) = p^a (1-p)^(1-a)
#  |         |
# a_1       a_2                             a_1     a_2
#  |         |                               |       |
# [ ]       [ ] p(s|a) = s^a (1-s)^(1-a)     v       v
#  |         |                              s_1 --> s_2
# s_1 -[ ]- s_2   
#  |    P(s2|s1) 
# [ ]
#  |    P(g|s) = 1.0
#  g    

sujetos = [1,2,3,4,5,6]

function pa_p(a,p=0.5)
    sum(a) == 0 ? p^2 : (sum(a) == 1 ? p*(1-p) : (1-p)^2)
end

function ps_a(s,a,e=1.5/2.1)
    a[mod(s-1,2)+1] == 1 ? e : 1-e
end

function grupo(s)
    div(s-1,2)
end

function ps_s(s1, s0)
    if grupo(s0) != grupo(s1)
        return 0.0
    elseif grupo(s0) == 0
        return 1/2
    elseif grupo(s0) == 1
        if s0 == 3
            return 1/2
        elseif s0 == s1
            return 1.0
        elseif s0 != s1 
            return 0.0
        end
    elseif grupo(s0) == 2
        if s0 == s1
            return 1.0
        else
            return 0.0
        end
    end
end

function ps2a1(s2,a1)
    res = 0.0
    for s1 in sujetos 
        res += pa_p(a1)*ps_a(s1,a1)*ps_s(s2,s1)
    end
    return res
end

function ps2_a1(s2,a1)
    tot = 0.0
    for s in sujetos 
        tot += ps2a1(s,a1)
    end
    return ps2a1(s2,a1)/tot
end

ps2_a1(1,[0,0])

function ps2a1a2(s2,a1,a2)
    ps2a1(s2,a1)*ps_a(s2,a2)
end

function ps2_a1a2(s2,a1,a2)
    tot = 0.0
    for s in sujetos    
        tot += ps2a1a2(s,a1,a2)
    end
    return ps2a1a2(s2,a1,a2)/tot
end

ps2_a1a2(4,[1,0],[0,1])

function pg_s(g,s)
    g == grupo(s) ? 1.0 : 0.0
end

function pga1a2(g,a1,a2)
    res = 0.0
    for s in sujetos
        res += ps2a1a2(s,a1,a2)*pg_s(g,s)
    end
    return res
end

function pg_a1a2(g,a1,a2)
    tot = 0.0
    for G in 0:2
        tot += pga1a2(G,a1,a2)
    end
    return pga1a2(g,a1,a2)/tot
end

pg_a1a2(0,[0,1],[1,0])
