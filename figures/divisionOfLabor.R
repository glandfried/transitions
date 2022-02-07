oldpar <- par(no.readonly = TRUE)
oldwd <- getwd()
this.dir <- dirname(parent.frame(2)$ofile)
nombre.R <-  sys.frame(1)$ofile
require(tools)
nombre <- print(file_path_sans_ext(nombre.R))
nombre <- "divisionOfLabor"
pdf(paste0("pdf/",nombre,".pdf"), width = 5, height = 5  )
setwd(this.dir)
###############################



base_2 <- function(n=2, pad=4){
    res = c()
    while (n > 1){
        res = c(res, n%%2 )
        n = n %/% 2
    }
    res = c(res, n )
    res = c(res, rep(0,pad-length(res)))
    return(res)
}

base_10 <- function(xs){
    res = 0
    for (i in seq(length(xs))){
        res = res + xs[i]*2^(i-1)
    }
    return(res)
}

p_as <- function(as, p){
    res = 1
    for (a in as){
        res = res * p^(a) * (1-p)^(1-a)
    }
    return(res)
}

p_as(c(1,1),p=0.75)

g_C <- function(ess, p = c(0.9,0.9)){
    res = 1
    N = dim(ess)[1]
    A = base_10(rep(1,N))
    D = dim(ess)[2]/2
    for (d in seq(0,D-1)){#d=1
        for (a in seq(0,A)){#a=0
            as = base_2(a,pad=N)
            common = 0
            for (i in seq(N)){#i=1
                common = common + (1/N)*(ess[i,(d*2)+1]^as[i]  * ess[i,(d*2)+2]^(1-as[i]) )
            }
            res = res * common^(p_as(as, p[d+1])/D)
        }
    }
    return(res)
}


distintas_estrategias <- function(a,b){
    res <- matrix(0,nrow=2,ncol=2)
    res[1,] <- c(a, 1-a)
    res[2,] <- c(b, 1-b)
    return(res)
} 

es = seq(0.01,0.99,by=0.01)
fitness = matrix(0, nrow = length(es), ncol=length(es))

for (i in seq(length(es))){
    for (j in seq(length(es))){
        fitness[i,j] = g_C(distintas_estrategias(es[i],es[j]), p=c(0.6) )
    }
}

image(fitness,axes = F,ann = F)
contour(fitness, add=T, levels=seq(0.1,0.6, by=0.02),  labcex=1)
#fitness[which.max(fitness)%/%dim(fitness)[1]+1, which.max(fitness)%%dim(fitness)[2]] == max(fitness)
points(es[which.max(fitness)%%dim(fitness)[2]],es[which.max(fitness)%/%dim(fitness)[1]+1], pch=19, cex=0.5)
axis(side=2, labels=NA,cex.axis=0.6,tck=0.015)
axis(side=1, labels=NA,cex.axis=0.6,tck=0.015)
axis(lwd=0,side=1, las=0,cex.axis=1.33,line=-0.45)
axis(lwd=0,side=2,cex.axis=1.33,line=-0.45)

mtext(text = "Estrategy 1",side =2 ,line=2,cex=1.66)
mtext(text = "Estrategy 2",side =1 ,line=2,cex=1.66)

dif <- seq(-0.5,0.5,by=0.02)
fitness_monopolio = matrix(NA, nrow = length(es), ncol= length(dif))
for (i in seq(length(es))){#i=1
    for (j in seq(length(dif))){#j=1
        ess = distintas_estrategias(es[i],es[i])
        ess[1, ] = ess[1,]+dif[j]
        ess[2,] = ess[2,]-dif[j]
        if (sum(ess<0) ==0){
            fitness_monopolio[i,j] = g_C(ess, p=c(0.6))
        }
    }
}

image(es, dif*2, fitness_monopolio,axes = F,ann = F, xlim=c(-0.01, 1.01))
contour(es, dif*2, fitness_monopolio, add=T, levels=seq(0.1,0.6, by=0.02),  labcex=1)
#fitness[which.max(fitness)%/%dim(fitness)[1]+1, which.max(fitness)%%dim(fitness)[2]] == max(fitness)
segments(x0 = 0.69, x1 = 0.69, y0 = 0.31*2, y1=-0.31*2, lwd=1.5, lty=3)
axis(side=2, labels=NA,cex.axis=0.6,tck=0.015)
axis(side=1, labels=NA,cex.axis=0.6,tck=0.015)
axis(lwd=0,side=1, las=0,cex.axis=1.33,line=-0.45)
axis(lwd=0,side=2,cex.axis=1.33,line=-0.45)

mtext(text = "Estrategy",side =1 ,line=2,cex=1.66)
mtext(text = "Division of labor",side =2 ,line=2,cex=1.66)


if (F){
    par(mar=c(3.75,3.75,0.25,0.25))


    GENERALIST <- matrix(0,nrow=4,ncol=4)
    GENERALIST[2,] <- c(0.45, 0.05, 0.45, 0.05)
    GENERALIST[3,] <- c(0.45, 0.05, 0.45, 0.05)
    GENERALIST[4,] <- c(0.45, 0.05, 0.45, 0.05)
    GENERALIST[1,] <- c(0.45, 0.05, 0.45, 0.05)
    g_C(GENERALIST,p)

    random_optimizer <- function(ess = GENERALIST , p = c(0.9, 0.9), alpha = 0.01, max_ITER = 1000){
        
        max_r = g_C(ess, p)
        iterations = 0
        step = Inf
        while ((iterations < max_ITER)){
            iterations = iterations + 1
            
            #row = sample(1:4,1)#row=4
            #cols = sample(1:4,2)#cols=c(1,3)
            
            desde = sample(which(ess > alpha + alpha/2),1)
            hacia = sample(which(ess < 1 - alpha - alpha/2),2)
            if (desde == hacia[1]){
                hacia = hacia[2]
            }else{
                hacia = hacia[1]
            }
            
            new_ess = ess
            new_ess[desde] = new_ess[desde] - alpha
            new_ess[hacia] = new_ess[hacia] + alpha
                        
            if (g_C(new_ess,p) - max_r > 1e-16){
                step =  g_C(new_ess,p) - max_r
                max_r = g_C(new_ess,p)
                ess = new_ess
                iterations = 0
            }
            
            

            
        }
        return(new_ess)
    }

    opt = random_optimizer(max_ITER=1000, alpha=0.01)
    g_C(opt,p=c(0.9, 0.9))

    E99_4346857 <- matrix(0,nrow=4,ncol=4)
    E99_4346857[1,]  <- c(0.49, 0.01, 0.49, 0.01)
    E99_4346857[2,]  <- c(0.49, 0.01, 0.49, 0.01)
    E99_4346857[3,]  <- c(0.49, 0.01, 0.49, 0.01)
    E99_4346857[4,]  <- c(0.49, 0.01, 0.49, 0.01)
    g_C(E99_4346857,c(0.9,0.9))

    opt = random_optimizer(p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    opt = random_optimizer(opt,p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    g_C(opt,p=c(0.9, 0.9))

    E65_2595482 <- matrix(0,nrow=4,ncol=4)
    E65_2595482[1,]  <- c(0.42, 0.07, 0.35, 0.15)
    E65_2595482[2,]  <- c(0.42, 0.07, 0.35, 0.15)
    E65_2595482[3,]  <- c(0.42, 0.07, 0.34, 0.15)
    E65_2595482[4,]  <- c(0.43, 0.08, 0.34, 0.15)
    g_C(E65_2595482,c(0.6,0.55))

    opt = random_optimizer(E65_2595482,p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    g_C(opt,c(0.6,0.55)) -g_C(E65_2595482,c(0.6,0.55))


    raro <- matrix(0,nrow=4,ncol=4)
    raro [1,]  <- c(1, 1, 0, 0)
    raro [2,]  <- c(1, 0, 1, 0)
    raro [3,]  <- c(0, 0, 1, 1)
    raro [4,]  <- c(0, 1, 0, 1)
    g_C(raro/2,c(0.6,0.55))

    opt = random_optimizer(raro/2,p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    g_C(opt,c(0.6,0.55)) -g_C(E65_2595482,c(0.6,0.55))
    opt = random_optimizer(opt,p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    g_C(opt,c(0.6,0.55)) -g_C(E65_2595482,c(0.6,0.55))

    c(sum(opt[,1]),sum(opt[,2]),sum(opt[,3]),sum(opt[,4]))/4


    raro2 <- matrix(0,nrow=4,ncol=4)
    raro2 [1,]  <- c(1, 1, 0, 0)
    raro2 [2,]  <- c(0, 0, 1, 1)
    raro2 [3,]  <- c(0, 0, 1, 1)
    raro2 [4,]  <- c(1, 1, 0, 0)
    g_C(raro2/2,c(0.6,0.55))
    opt = random_optimizer(raro2/2,p=c(0.6,0.55),max_ITER=1000, alpha=0.01)
    g_C(opt,c(0.6,0.55)) -g_C(E65_2595482,c(0.6,0.55))

    c(sum(opt[,1]),sum(opt[,2]),sum(opt[,3]),sum(opt[,4]))/4

    c(sum(opt[1,]),sum(opt[2,]),sum(opt[3,]),sum(opt[4,]))


    siempre_opt <- matrix(0,nrow=4,ncol=4)
    siempre_opt[1,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)
    siempre_opt[2,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)
    siempre_opt[3,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)
    siempre_opt[4,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)
    g_C(siempre_opt,c(0.6,0.55))



    contour_opt  <- matrix(0,nrow=4,ncol=4)
    contour_opt [1,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)-0.0775
    contour_opt [2,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)+0.2325
    contour_opt [3,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)-0.0775
    contour_opt [4,]  <- c(0.4250, 0.0775, 0.3450, 0.1525)-0.0775
    g_C(siempre_opt,c(0.6,0.55)) - g_C(contour_opt ,c(0.6,0.55))



    contour_opt <- matrix(0,nrow=4,ncol=4)
    contour_opt[1,]  <- c(0.55, 0.20, 0.20, 0.01)
    contour_opt[2,]  <- c(0.39, 0.05, 0.26, 0.06)
    contour_opt[3,]  <- c(0.36, 0.01, 0.57, 0.38)
    contour_opt[4,]  <- c(0.40, 0.05, 0.35, 0.16)
    g_C(siempre_opt,c(0.6,0.55)) - g_C(contour_opt,c(0.6,0.55))

    dos <- matrix(0,nrow=2,ncol=2)
    dos[1,]  <- c(0.1,0.9)
    dos[2,]  <- c(0.9,0.1)
    opt = random_optimizer(dos,p=c(0.6),max_ITER=1000, alpha=0.01)
    opt = random_optimizer(opt,p=c(0.6),max_ITER=1000, alpha=0.01)

    dos <- matrix(0,nrow=2,ncol=2)
    dos[1,]  <- c(0.9,0.1)
    dos[2,]  <- c(0.1,0.9)

    opt[,1]/opt[,2]
    sum(opt[,1])/2
    sum(opt[,2])/2

    dos <- matrix(0,nrow=2,ncol=2)
    dos[1,]  <- c(0.69,0.31)
    dos[2,]  <- c(0.69,0.31)
    g_C(dos,p=(0.6))

    contour_dos <- matrix(0,nrow=2,ncol=2)
    contour_dos[1,]  <- c(0.38,0.00)
    contour_dos[2,]  <- c(1,0.62)
    g_C(dos,p=(0.6)) - g_C(contour_dos,p=(0.6))


    (0.60+0.22)
    0.78/(0.78+0.40)

}


#######################################
# end 
dev.off()
system(paste("pdfcrop -m '0 0 0 0'",paste0("pdf/",nombre,".pdf") ,paste0("pdf/",nombre,".pdf")))
setwd(oldwd)
par(oldpar, new=F)
#########################################
