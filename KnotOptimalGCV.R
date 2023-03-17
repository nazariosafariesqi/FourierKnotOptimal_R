install.packages("MASS")
library(MASS)
setwd("C:/Users/Nazario/Desktop")
getwd()
reg=read.csv("kemiskinan.csv", header = T)
reg
y= as.matrix(reg[,1])
y
x1=as.matrix(reg[,2])
x1
x2=as.matrix(reg[,3])
x2
x=data.frame(y,x1,x2)
x
x=as.matrix(x)
x

MPL=function(x,eps=1e-0009)
{
  x=as.matrix(x)
  xsvd=svd(x)
  diago=xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus=as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u)/diago)
  }
  else
  {
    xplus=xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}
MPL(reg,eps=1e-0009)
deretfourier=function(y,x1,x2,K)
{
  n=length(y)
  a=(4*(K+1))+1
  H=matrix(0,n,a)
  hasil=matrix(0,K,2)
  for(k in 1:K)
  {
    for(i in 1:n)
    {
      for(j in 1:k)
      {
        H[i,1]=1
        H[i,2]=x1[i]
        H[i,2+j]=cos(j*x1[i])
        H[i,3+k]=x2[i]
        H[i,3+j+k]=cos(j*x2[i])
      }
    }
    I=diag(1,n,n)
    A=H%*%MPL(t(H)%*%H)%*%t(H)
    Atot=A
    ytop=Atot%*%y
    df=sum(diag(Atot))
    W=(1/n)*I
    atas=t(y-ytop)%*%W%*%(y-ytop)
    bawah=((1-df)/n)^2
    GCV=atas/bawah
    hasil[k,1]=k
    hasil[k,2]=GCV
  }
  print(hasil)
  
  GCV2=min(hasil[,2])
  s=1
  repeat{
    if(hasil[s,2]==GCV2)
    {
      k.optimal=hasil[s,1]
      GCVOpt=GCV2
      break
    }
    else s=s+1
  }
  cat("Nilai K optimal adalah \t", k.optimal, "\n")
  print(H)
}

