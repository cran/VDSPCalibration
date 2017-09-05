samplesize<- function (x0, d0,cutpts=c(7.5,42.5,57.5,72.5,200),CVx, CVy){

    x=seq(cutpts[1],cutpts[2],length=1000)
    i = 2
    while(i < length(cutpts)){
      x = c(x, seq(cutpts[i],cutpts[i+1], length=1000))
      i = i + 1
    }

    mu1 = mean(x)
    mu2 = mean(x^2)/(1 + CVx^2)
    mu3 = mean(x^3)/(1 + 3 * CVx^2)
    mu4 = mean(x^4)/(1 + 6 * CVx^2 + 3 * CVx^4)
    v11 = (CVx^2 + CVy^2) * mu2
    v12 = (CVy^2 + (CVx^2 - CVx^4)/(1 + CVx^2)) * mu3
    v13 = (CVx^2 + (CVy^2 - CVy^4)/(1 + CVy^2)) * mu3
    v22 = (CVy^2 * (1 + CVx^2) + CVx^2 * (1 + CVx^4)/(1 + CVx^2)^2) *mu4
    v23 = ((CVx^2 - CVx^4)/(1 + CVx^2) + (CVy^2 - CVy^4)/(1 +
        CVy^2) + CVx^2 * CVy^2/(1 + CVx^2)/(1 + CVy^2)) * mu4
    v33 = (CVx^2 * (1 + CVy^2) + CVy^2 * (1 + CVy^4)/(1 + CVy^2)^2) *mu4
    A = cbind(c(1, mu1), c(mu1, mu2), c(mu1, mu2))
    B = cbind(c(v11, v12, v13), c(v12, v22, v23), c(v13, v23, v33))
    sigma = solve(A %*% solve(B) %*% t(A))
    x0 = c(1, x0)
    v = t(x0) %*% sigma %*% x0
    n = (2 * 1.96)^2/d0^2 * v
    size = ceiling(n)[1,1]
    return(size = size)
}

samplefun=function(x, index, n0)
{
id=order(x)
x=x[id]
index=index[id]
m=length(x)
x.target=seq(x[1], x[m], length=n0)

x.select=rep(0, n0)
index.select=rep(0,n0)

x.candidate=x
index.candidate=index

for(i in 1:n0)
    {distance=abs(x.candidate-x.target[i])
     m0=length(x.candidate)
     id0=(1:m0)[distance==min(distance)][1]

     index.select[i]=index.candidate[id0]
     x.select[i]=x.candidate[id0]

     index.candidate=index.candidate[-id0]
     x.candidate=x.candidate[-id0]
     }
return(list(index=index.select, x=x.select))
}


sampletot=function(x, index, n0, K)
  {m0=round(n0/K)
   cut=quantile(x, seq(0, 1, length=K+1))
   id.result=NULL
   id.x=NULL
   for(k in 1:K)
    {xsub=x[x>=cut[k] & x<cut[k+1]]
     idsub=index[x>=cut[k] & x<cut[k+1]]
     fit=samplefun(xsub, idsub, m0)
     id.result=c(id.result, fit$index)
     id.x=c(id.x, fit$x)
     }
   return(list(x=id.x, index=id.result))
   }

calfun<- function (x, y, CVx, CVy = CVx, lambda0 = 1)
{
    n = length(x)
    A1 = cbind(c(1, mean(x)), c(mean(x), mean(x^2/(1 + CVx^2))))
    beta1 = solve(A1) %*% c(mean(y), mean(x * y))
    s1 = rbind(as.vector(y - beta1[1] - beta1[2] * x), x * as.vector(y - 
                             beta1[1] - beta1[2] * x/(1 + CVx^2)))
    B1 = s1 %*% t(s1)/n
    v1 = solve(A1) %*% B1 %*% solve(A1)/n
    sigma1 = sqrt(diag(v1))
    A2 = cbind(c(1, mean(y)), c(mean(x), mean(x * y)))
    beta2 = solve(A2) %*% c(mean(y), mean(y^2/(1 + CVy^2)))
    s2 = rbind(as.vector(y - beta2[1] - beta2[2] * x), y * as.vector(y/(1 + 
                         CVy^2) - beta2[1] - beta2[2] * x))
    B2 = s2 %*% t(s2)/n
    v2 = solve(A2) %*% B2 %*% solve(A2)/n
    sigma2 = sqrt(diag(v2))
    s3 = rbind(as.vector(y - beta2[1] - beta2[2] * x), x * as.vector(y - 
                             beta2[1] - beta2[2] * x/(1 + CVx^2)), 
                             y * as.vector(y/(1 +  CVy^2) - beta2[1] - 
                             beta2[2] * x))
    sigma3 = solve(s3 %*% t(s3))
    objq = function(beta, sigma) {
        score = rbind(as.vector(y - beta[1] - beta[2] * x), x *
            as.vector(y - beta[1] - beta[2] * x/(1 + CVx^2)),
            y * as.vector(y/(1 + CVy^2) - beta[1] - beta[2] *
                x))
        score = apply(score, 1, sum)
        result = t(score) %*% sigma %*% score
        return(result)
    }
    beta3 = optim(beta2, objq, sigma = sigma3)$par
    s3 = rbind(as.vector(y - beta3[1] - beta3[2] * x), x * as.vector(y -
                             beta3[1] - beta3[2] * x/(1 + CVx^2)), 
                             y * as.vector(y/(1 + CVy^2) - beta3[1] - 
                             beta3[2] * x))
    sigma3 = solve(s3 %*% t(s3))
    beta3 = optim(beta2, objq, sigma = sigma3)$par
    s3 = rbind(as.vector(y - beta3[1] - beta3[2] * x), x * as.vector(y - 
                             beta3[1] - beta3[2] * x/(1 + CVx^2)), 
                             y * as.vector(y/(1 + CVy^2) - beta3[1] - 
                             beta3[2] * x))
    B3 = s3 %*% t(s3)/n
    A3 = cbind(c(1, mean(x)), c(mean(x), mean(x^2/(1 + CVx^2))),
        c(mean(y), mean(x * y)))
    v3 = solve(A3 %*% solve(B3) %*% t(A3))/n
    sigma3 = sqrt(diag(v3))
    beta4 = rep(0, 2)
    a = cov(x, y) * (lambda0 * var(x) + mean(x)^2)
    b = (cov(x, y)^2 - mean(x)^2 * var(y)) - lambda0 * (cov(x, y)^2 - 
         var(x) * mean(y)^2)
    c = -cov(x, y) * (var(y) + lambda0 * mean(y)^2)
    beta4[2] = (-b + sqrt(b^2 - 4 * a * c))/2/a
    beta4[1] = mean(y) - beta4[2] * mean(x)
    rx2 = beta4[2] * mean(x^2)/(mean(x * y) - beta4[1] * mean(x))
    ry2 = mean(y^2)/(beta4[1] * mean(y) + beta4[2] * mean(x * y))
    A4 = cbind(c(1, rx2 * mean(x), ry2 * mean(y)), c(mean(x),
        mean(x^2), ry2 * mean(x * y)), c(0, beta4[1] * mean(x) -
        mean(x * y), lambda0 * beta4[1] * mean(y) + lambda0 *
        beta4[2] * mean(x * y)))
    s4 = rbind(as.vector(y - beta4[1] - beta4[2] * x), x * as.vector(y *rx2 - 
                             beta4[1] * rx2 - beta4[2] * x), y * as.vector(y -  
                             beta4[1] * ry2 - beta4[2] * x * ry2))
    B4 = s4 %*% t(s4)/n
    v4 = (solve(A4) %*% B4 %*% t(solve(A4))/n)[1:2, 1:2]
    sigma4 = sqrt(diag(v4))
    betax = beta1
    betay = beta2
    betaxy = beta3
    betaratio = beta4
    sex = sigma1
    sey = sigma2
    sexy = sigma3
    seratio = sigma4
    result = matrix(0, 4, 8)
    result[1, ] = c(betax, betay, betaxy, betaratio)
    result[2, ] = c(sex, sey, sexy, seratio)
    result[3, ] = c(betax - 1.96 * sex, betay - 1.96 * sey, 
                    betaxy -1.96 * sexy, betaratio - 1.96 * seratio)
    result[4, ] = c(betax + 1.96 * sex, betay + 1.96 * sey, 
                    betaxy + 1.96 * sexy, betaratio + 1.96 * seratio)
    colnames(result) = c("intercept-CVx only", "slope-CVx only",
        "intercept-CVy only", "slope-CVy only", "intercept-CVx and CVy",
        "slope-CVx and CVy", "intercept-CVy/CVx only", "slope-CVy/CVx")
    rownames(result) = c("coef", "se", "lower CI", "upper CI")
    return(result)
}


chngpt = function(x,y,start=quantile(x,probs=0.10,na.rm="TRUE"), 
                  finish=quantile(x,probs=0.90,na.rm="TRUE"), NbrSteps=500){

 yfit = lm( y ~ x)
 rms_min = summary(yfit)$sigma
 xval = NA
 coefficients = summary(yfit)$coefficients
 yfitted = yfit$fitted.values

 d = (finish - start)/NbrSteps

 tau = start
 while (tau < finish){
    x1 = x
    x2 = (x - tau)*(x >= tau)
    yfit = lm(y~x + x2)
    rms = summary(yfit)$sigma
    if (rms < rms_min) {
       rms_min = rms
       xval = tau
       coefficients = summary(yfit)$coefficients
       yfitted = yfit$fitted.values
       }
    tau = tau + d
 }
return(list(x=x, y=y, yfitted=yfitted, chngpt=xval, coefficients=coefficients))
}

