
####################################### function ####################################
########
############熵权法
########

CWfunc <- function(data,type){
  ##data:收集到的数据矩阵，行名是指标，列名是年份
  ##type:正指标为1；负指标为0
  ##
  ##
  sum(is.na(data))
  ###归一化处理
  min.max.norm <- function(x){
    (x-min(x))/(max(x)-min(x))
  }
  
  max.min.norm <- function(x){
    (max(x)-x)/(max(x)-min(x))
  }
  
  data_1 = data[lab==1,]
  data_1 = apply(data_1,2,min.max.norm)  #正向
  data_0 = data[lab==0,]
  data_0 = apply(data_0,2,max.min.norm)   #负向
  
  data <- cbind(data_1,data_0)
  
  ###求出所有样本对指标Xj的贡献总量  
  func1 <- function(data)
  {
    x <- c(data)
    for(i in 1:length(data))
      x[i] = data[i]/sum(data[])
    return(x)
  }
  Con_total <- apply(data,2,func1)
  
  ###将上步生成的矩阵每个元素变成每个元素与该ln（元素）的积并计算信息熵 
  fun2 <- function(data)
  {
    x <- c(data)
    for(i in 1:length(data)){
      if(data[i] == 0){
        x[i] = 0
      }else{
        x[i] = data[i] * log(data[i])
      }
    }
    return(x)
  }
  IE <- apply(Con_total,2,fun2)
  
  k <- 1/log(length(IE[,1]))
  d <- -k * colSums(IE)
  ###计算冗余度
  d <- 1-d
  ###计算权重
  w <- d/sum(d)
  return(w)
}



#############
####灰熵模型函数部分
###############


###############
###error correction model 
##############

ECMfunc <- function(data,k,C){
  #将特定数据进行区间化
  #data: 输入的数据框结构的样本,行为年份,列为指标;
  #k: 为误差参数,随实际情况进行调节,一般低于0.01%;
  #C: 为均匀取样参数,一般大于10000;
  ECM_data=list();
  ECM_data.temp=list();
  temp=NULL;
  for (i in 1:len_i) {
    for (j in 1:len_j) {
      range_low=data[i,j]*(1-0.01*k)
      range_high=data[i,j]*(1+0.01*k)
      temp = runif(C,min=range_low,max=range_high)
      ECM_data.temp[[j]]=temp
    }
    ECM_data[[i]]=ECM_data.temp
  }
  return(ECM_data)
}


###########################################
###############
###error correction border 
##############

ECBfunc <- function(data,k){
  #获得特定数据进行区间化的数据边界
  #data: 输入的数据框结构的样本,行为年份,列为指标;
  #k: 为误差参数,随实际情况进行调节,一般低于0.01%;
  ECB_border=list();
  ECB_border.temp=list();
  temp=NULL;
  for (i in 1:len_i) {
    for (j in 1:len_j) {
      range_low=data[i,j]*(1-0.01*k)
      range_high=data[i,j]*(1+0.01*k)
      temp = c(min=range_low,max=range_high)
      ECB_border.temp[[j]]=temp
    }
    ECB_border[[i]]=ECB_border.temp
  }
  return(ECB_border)
}

###############
###Remove dimensional effects
##############

RDEminfunc <-function(data,ECB_data){
  #计算上边界去量纲分母
  #ECB_data:为经过误差校正(ECBfunc)得到的数据边界值,为list结构;
  sum_min=NULL
  for (i in 1:len_i) {
    for (j in 1:len_j) {
      sum_temp_min=sqrt(sum((ECB_data[[i]][[j]]["min"])^2))
      }
    sum_min=c(sum_min,sum_temp_min)
  }
  return(sum_min)
}

RDEmaxfunc <-function(data,ECB_data){
  #计算下边界去量纲分母
  #ECB_data:为经过误差校正(ECBfunc)得到的数据边界值,为list结构;
  sum_max=NULL
  for (i in 1:len_i) {
    for (j in 1:len_j) {
      sum_temp_max=sqrt(sum((ECB_data[[i]][[j]]["max"])^2))
    }
    sum_max=c(sum_max,sum_temp_max)
  }
  return(sum_max)
}

RDEfunc <- function(ECM_data,data_min,data_max,C){
  #去除量纲影响
  #ECM_data:为经过误差校正(ECMfunc)得到的数据,为list结构;
  #data_min:为经过函数(RDEminfunc)得到的数据;
  #data_max:为经过函数(RDEmaxfunc)得到的数据;
  #C: 为均匀取样参数,一般大于10000;
  len_k = C
  RDE_data_temp1=list();
  RDE_data_temp2=list();
  RDE_data = list();
  
  for (i in 1:len_i) {
    for (j in 1:len_j ) {
      for (k in 1:len_k) {
          RDE_data_temp_max = NULL;
          RDE_data_temp_min = NULL;
          RDE_data_temp_max = ECM_data[[i]][[j]][k]/data_min[i];  #
          RDE_data_temp_min = ECM_data[[i]][[j]][k]/data_max[i];  #
          RDE_data_temp1[[k]]= c(RDE_data_temp_min,RDE_data_temp_max)
      }
      RDE_data_temp2[[j]]=RDE_data_temp1
    }
    RDE_data[[i]]=RDE_data_temp2
  }
  return(RDE_data)
}


############计算参考序列
SSEQ <- function(RDE_data,C){
  ######此函数计算参考数列
  ######RDE_data:为经过函数(RDEfunc)去除量纲影响后得到的数据;
  ######C: 为均匀取样参数,一般大于10000;
  len_k = C;
  temp_up = NULL;
  temp_low = NULL;
  Range = list();
  for (i in 1:len_i) {
    temp_max = data.frame();
    temp_min = data.frame();
    for (j in 1:len_j ) {
      for (k in 1:len_k) {
        temp_max[j,k] = RDE_data[[i]][[j]][[k]][1]
        temp_min[j,k] = RDE_data[[i]][[j]][[k]][2]
      }
    }
    temp_up = sapply(temp_max,function(x){max(x)})       #参考序列上界
    temp_low = sapply(temp_min,function(x){max(x)})      #参考序列下界
    Range[[i]] = cbind(temp_low,temp_up)
  }
  return(Range)
}


####################计算Euclidean距离
EDfunc <- function(RDE_data,Range,C){
  #此函数是计算Euclidean距离
  ####RDE_data:为经过函数(RDEfunc)去除量纲影响后得到的数据;
  #Range是经过SSEQ函数得到的list型参考序列
  #C: 为均匀取样参数,一般大于10000;
  len_k = C;
  ED=list();
  ED_temp = list();
  ED_temp1 = list();
  for (i in 1:len_i) {
    for (j in 1:len_j ) {
      ED_temp1 = NULL;
      for (k in 1:len_k) {
        x = NULL;
        x = sqrt((RDE_data[[i]][[j]][[k]][1]-Range[[i]][k,1])^2
                                    +(RDE_data[[i]][[j]][[k]][2]-Range[[i]][k,2])^2/2);  #Euclidean距离公式
        ED_temp1 = c(ED_temp1,x)
      }
      ED_temp[[j]]=ED_temp1
    }
    ED[[i]]=ED_temp
  }
  return(ED)
}





########计算区间数属性值与理想区间数属性值的关联系数
CCfunc <- function(ED_data,rou,C=10000){
  ##此函数计算区间数属性值与理想区间数属性值的关联系数
  ###ED_data：经过函数EDfunc计算出的Euclidean距离数据;
  ##rou：为设定值,范围为[0,1],此处去0.5;
  ##C: 为均匀取样参数,一般大于10000;
  len_k = C;
  CC=list();
  CC_temp = list();
  for (i in 1:len_i) {
    for (j in 1:len_j ) {
      CC_temp1 = NULL
      ED_min=min(sapply(ED_data[[i]],function(x){min(x)}))
      ED_max=max(sapply(ED_data[[i]],function(x){max(x)}))
      for (k in 1:len_k) {
        x = NULL;
        x = (ED_min+rou*ED_max)/(ED_data[[i]][[j]][[k]]+rou*ED_max);  #关联系数公式
        CC_temp1 = c(CC_temp1,x)
      }
      CC_temp[[j]]=CC_temp1
    }
    CC[[i]]=CC_temp
  }
  return(CC)
}



#############结构转换

GRDTfunc <- function(CC_data){
  ##为了让运行速度加快，
  ##将CC_data的结构内部结构转化为matrix
  ##CC_data：通过函数CCfunc计算出的区间数属性值与理想区间数属性值的关联系数
  GRDT=list();
  for (i in 1:len_i) {
    temp = NULL;
    for (j in 1:len_j){
     x = CC_data[[i]][[j]]
     temp=cbind(temp,x)       #结构转换
    }
    GRDT[[i]]=temp
  }
  return(GRDT)
}


################计算灰色关联度
GRDfunc <- function(GRDT_data,data){
  ##函数功能为计算灰色关联度
  ##GRDT_data:结构转换过的距离数据,结构为list
  ##data:原始数据
  GRD = NULL;
  for (i in 1:len_i) {
    temp = data.frame();
    temp = apply(GRDT_data[[i]],1,function(x){mean(x)})     #计算灰色关联度
    GRD=cbind(GRD,temp)
  }
  colnames(GRD) = rownames(data)
  return(GRD)
}


#############计算灰熵权重
GEWfunc <- function(GRD_data){
  ##计算灰熵的最终权重
  ##GRD_data为通过函数GRDfunc得到的灰色关联度，结构为matrix;
  GEW=list()
  GE=NULL;
  GE=apply(GRD_data,2,function(x){-sum(x*log(x))})       #计算灰熵权重
  W=GE/sum(GE)
  GEW=cbind(weight=W,GreyEntropy=GE)
  return(GEW)
}



###########==========
##########===========
####======================= main code ==========================###
##################

data <- read.csv("the factor of hefei.csv ",header = T, row.names=1)
data <- t(data)
len_i = dim(data)[1];
len_j = dim(data)[2];
######数据偏差纠正
ECM_data <- ECMfunc(data,0.01,10000)

##################去量纲化
ECB_data <- ECBfunc(data,0.01)
data_min = RDEminfunc(data,ECB_data)
data_max = RDEmaxfunc(data,ECB_data)
RDE_data = RDEfunc(ECM_data,data_min,data_max,10000)

#########计算参考序列
Range <- SSEQ(RDE_data,10000)

####################计算Euclidean距离
ED_data <- EDfunc(RDE_data,Range,10000)

####rou对照表
rou_low <- rbind(seq(0.1,0.5,0.1),c(0.091,0.167,0.231,0.286,0.333))
rou_up <- rbind(seq(0.6,1,0.1),c(0.375,0.412,0.444,0.474,0.5))
rownames(rou_low)=c("index","number")
rownames(rou_up)=c("index","number")


##################计算关联度
CC_data <- CCfunc(ED_data,rou = 0.5,C=10000)
GRDT_data <- GRDTfunc(CC_data)
GRD_data = GRDfunc(GRDT_data,data)

############计算灰熵权重
weight <- GEWfunc(GRD_data)
##############保存数据
save(weight,file="weight.Rdata")
#load(file="weight.Rdata")
#显示weight
#weight


############
####=========================熵权法主程序===============================###


data <- read.csv("the factor of hefei.csv ",header = T, row.names=1)
pre <- read.csv("pre.csv",header = T, row.names=1)
lab <- c(1,1,1,1,1)
W=CWfunc(data,lab)
##############保存数据
save(W,file="weight.Rdata")
#load(file="W.Rdata")
#显示weight
W

##############
###============================灰度预测======================================##
### load packages to be used, if not installed, please use 
install.packages("forecast")
install.packages("zoo")

library(forecast)
library(zoo)

###加入权重进行分析

dataWfunc1 <- function(data,weight){
  data_W = NULL;
  for (i in 1:len_i) {
    data_temp=NULL;
    for (j in 1:len_j) {
      X=NULL;
      X=data[j,i]*weight[i,1]
      data_temp=rbind(data_temp,X)
    }
   data_W=cbind(data_W,data_temp)
  }
  data_W=as.data.frame(data_W)
  colnames(data_W)=colnames(data)
  rownames(data_W)=rownames(data)
  return(data_W)
}

dataWfunc2 <- function(data,weight){
  data_W = NULL;
  for (i in 1:len_i) {
    data_temp=NULL;
    for (j in 1:len_j) {
      X=NULL;
      X=data[j,i]*weight[i]
      data_temp=rbind(data_temp,X)
    }
    data_W=cbind(data_W,data_temp)
  }
  data_W=as.data.frame(data_W)
  colnames(data_W)=colnames(data)
  rownames(data_W)=rownames(data)
  return(data_W)
}

data_W = dataWfunc(data,weight)


View(state.x77)
cor(state.x77)
install.packages("car")
library(car)

data_re=cbind(data_W,pre)
data_re1=cbind(data,pre)
cor(data_re)
scatterplotMatrix(data_re,main = "Scatter Plot Matrix")
fit <- lm(data_re$垃圾年产量~data_re$城镇人口+data_re$人均生活消费支出+data_re$平均住宅面积
          +data_re$合肥生产总值+data_re$住宅使用总面积,data=data_re)
summary(fit)

fit1 <- lm(data_re1$垃圾年产量~data_re1$城镇人口+data_re1$人均生活消费支出+data_re1$平均住宅面积
          +data_re1$合肥生产总值+data_re1$住宅使用总面积,data=data_re1)
summary(fit1)

step(fit1)
par(mfrow=c(2,2))
plot(fit1)


par(mfrow=c(1,1)) #设置画面布局
par(mfrow=c(1,2)) #设置画面布局
# 预测计算
dfp<-predict(fit1,type="response")
# 打印预测时
head(dfp,10)
# 合并数据
Data$垃圾年产量
X = as.numeric(names(dfp)) 
# 画图
plot(X,dfp,type = "o",col="blue",lty=1)
lines(Data$垃圾年产量,type = "o",col="black",lty=2)


plot(X,Data$垃圾年产量,type = "o",col="red",lty=3)
lines(dfp,type = "o",col="blue",lty=1)


###########
###误差分析

#plot(Test_fisher_ACC[4:6],type = "o",lty=1,xlab = "不同参数",ylab = "Error_Rate")
#lines(Test_stouffer_ACC[4:6],type = "o",col="yellow",lty=2)

######=================
############
########
head(data_re)
Data=data_re
###数据概括性度量
Min=sapply(Data,min)      #最小值
Max=sapply(Data,max)      #最大值
Mean=sapply(Data,mean)    #均值
SD=sapply(Data,sd)        #方差
cbind(Min,Max,Mean,SD)


#pearson相关系数，保留两位小数
round(cor(Data,method = c("pearson")),2)

## Loading required package: 
install.packages("lars")
library(lars)
#加载adapt-lasso源代码
lasso.adapt.bic2<-function(x,y){
  require(lars)
  ok<-complete.cases(x,y)
  x<-x[ok,]                           
  y<-y[ok]                             
  m<-ncol(x)
  n<-nrow(x)
  x<-as.matrix(x)                      
  one <- rep(1, n)
  meanx <- drop(one %*% x)/n
  xc <- scale(x, meanx, FALSE)         
  normx <- sqrt(drop(one %*% (xc^2)))
  names(normx) <- NULL
  xs <- scale(xc, FALSE, normx)        
  out.ls=lm(y~xs)                      
  beta.ols=out.ls$coeff[2:(m+1)]       
  w=abs(beta.ols)                      
  xs=scale(xs,center=FALSE,scale=1/w) 
  object=lars(xs,y,type="lasso",normalize=FALSE)
  sig2f=summary(out.ls)$sigma^2       
  bic2=log(n)*object$df+as.vector(object$RSS)/sig2f       
  step.bic2=which.min(bic2)            
  fit=predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
  coeff=predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
  coeff=coeff*w/normx                 
  st=sum(coeff !=0)                    
  mse=sum((y-fit)^2)/(n-st-1)          
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
  intercept=as.numeric(mean(y)-meanx%*%coeff)
  return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,intercept=intercept,object=object,
              bic2=bic2,step.bic2=step.bic2))
}


out1<-lasso.adapt.bic2(x=Data[,1:5],y=Data$垃圾年产量)


#adapt-lasso输出结果名称
names(out1)

#变量选择输出结果序号
out1$x.ind

#保留五位小数
round(out1$coeff,5)


#加载GM(1,1)源文件
gm11<-function(x0,t){ #x0为输入学列，t为预测个数
  x1<-cumsum(x0) #一次累加生成序列1-AG0序列
  b<-numeric(length(x0)-1)
  n<-length(x0)-1
  for(i in 1:n){ #生成x1的紧邻均值生成序列
    b[i]<--(x1[i]+x1[i+1])/2
    b} #得序列b，即为x1的紧邻均值生成序列
  D<-numeric(length(x0)-1)
  D[]<-1
  B<-cbind(b,D)
  BT<-t(B)#做逆矩阵
  M<-solve(BT%*%B)
  YN<-numeric(length(x0)-1)
  YN<-x0[2:length(x0)]
  alpha<-M%*%BT%*%YN  #模型的最小二乘估计参数列满足alpha
  alpha2<-matrix(alpha,ncol=1)
  a<-alpha2[1]
  u<-alpha2[2]
  cat("GM(1,1)参数估计值：",'\n',"发展系数-a=",-a,"  ","灰色作用量u=",u,'\n','\n') #利用最小二乘法求得参数估计值a,u
  y<-numeric(length(c(1:t)))
  y[1]<-x1[1]
  for(w in 1:(t-1)){  #将a,u的估计值代入时间响应序列函数计算x1拟合序列y
    y[w+1]<-(x1[1]-u/a)*exp(-a*w)+u/a 
  }
  ##cat("x(1)的模拟值：",'\n',y,'\n')
  xy<-numeric(length(y))
  xy[1]<-y[1]
  for(o in 2:t){ #运用后减运算还原得模型输入序列x0预测序列
    xy[o]<-y[o]-y[o-1] 
  } 
  cat("x(0)的模拟值：",'\n',xy,'\n','\n')                       
  
  #计算残差e
  e<-numeric(length(x0))
  for(l in 1:length(x0)){
    e[l]<-x0[l]-xy[l] #得残差
  }
  ##cat("残差：",'\n',e,'\n')
  #计算相对误差
  e2<-numeric(length(x0))
  for(s in 1:length(x0)){
    e2[s]<-(abs(e[s])/x0[s]) #得相对误差
  }
  ##cat("相对残差：",'\n',e2,'\n','\n')
  cat("残差平方和=",sum(e^2),'\n')
  cat("平均相对误差=",sum(e2)/(length(e2)-1)*100,"%",'\n')
  cat("相对精度=",(1-(sum(e2)/(length(e2)-1)))*100,"%",'\n','\n')
  
  #后验差比值检验
  avge<-mean(abs(e));esum<-sum((abs(e)-avge)^2);evar=esum/(length(e)-1);se=sqrt(evar)  #计算残差的方差se
  avgx0<-mean(x0);x0sum<-sum((x0-avgx0)^2);x0var=x0sum/(length(x0));sx=sqrt(x0var)  #计算原序列x0的方差sx
  cv<-se/sx  #得验差比值
  cat("后验差比值检验:",'\n',"C值=",cv,'\n')#对后验差比值进行检验，与一般标准进行比较判断预测结果好坏。
  if(cv < 0.35){     
    cat("C值<0.35, GM(1,1)预测精度等级为：好",'\n','\n')
  }else{
    if(cv<0.5){
      cat("C值属于[0.35,0.5), GM(1,1)模型预测精度等级为：合格",'\n','\n')
    }else{
      if(cv<0.65){
        cat("C值属于[0.5,0.65), GM(1,1)模型预测精度等级为：勉强合格",'\n','\n')
      }else{
        cat("C值>=0.65, GM(1,1)模型预测精度等级为：不合格",'\n','\n')
      }
    }
  }
  #画出输入序列x0的预测序列及x0的比较图像
  plot(xy,col='blue',type='b',pch=16,xlab='时间序列',ylab='值')
  points(x0,col='red',type='b',pch=4)
  legend('topleft',c('预测值','原始值'),pch=c(16,4),lty=l,col=c('blue','red'))
}

gm11(Data$城镇人口/10000,length(Data$城镇人口/10000)+6)
gm11(Data$人均生活消费支出,length(Data$人均生活消费支出)+2)
gm11(Data$合肥生产总值,length(Data$合肥生产总值)+2)
gm11(Data$平均住宅面积,length(Data$平均住宅面积)+2)
gm11(Data$住宅使用总面积,length(Data$住宅使用总面积)+2)
gm11(Data$垃圾年产量[1:10],length(Data$垃圾年产量[1:10])+3)

Data[21,1:6]
X <- Data
X[20,1:6] <- NA
X <- na.omit(X)
gm11(X$垃圾年产量,length(X$垃圾年产量)+3)

XX =gm11(data_re1$城镇人口,length(data_re1$垃圾年产量)+2)

X1=gm11(Data$城镇人口,length(Data$垃圾年产量)+2)


#读入数据
library(nnet)
#head(Data)
asData=scale(Data)
y=Data$垃圾年产量
nn<-nnet(y~.,asData[1:23,1:4],size=4,decay=0.0000001,maxit=10000,linout=T,trace=T)
predict<-predict(nn,asData[,1:6])
a=1:23
#画出序列预测值、真实值图像
plot(predict[1:23,1],col='red',type='b',pch=16,xlab='年份',ylab='垃圾年产量',xaxt="n")
points(Data[1:23,6],col='blue',type='b',pch=4)
legend('topleft',c('垃圾年产量预测值','垃圾年产量真实值'),pch=c(1,23),col=c('red','blue'))
axis(1,at=1:23,labels=a)
