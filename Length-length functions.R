
#Functions "r2" and "fn.get.TL" compute the r2 of the robust linear regression and ....

library(foreign)
library(MASS)


#r2 funtion:

r2 <- function(x)
{
  SSe <- sum((x$resid)^2)
  observed <- x$resid+x$fitted
  SSt <- sum((observed-mean(observed))^2)
  value <- 1-SSe/SSt
  return(value)
}


#fn.get.TL function:

fn.get.TL=function(DAT,SP,Show.Text,N,By.sex)
{
  test=subset(DAT,COMMON_NAME%in%SP & TL>0)
  for(i in 1:length(SP))
  {
    b=subset(test,COMMON_NAME==SP[i])
    if(Show.Text=="YES")
    {
      plot(b$TL,b$FL,main=SP[i],col="transparent")
      text(b$TL,b$FL,labels=b$YEAR,offset=0.5,cex=0.85)
    }else{
      plot(b$TL,b$FL,main=SP[i],col=1,pch=19)
    }
    if(nrow(b)>=N)
    {
      if(By.sex=="YES")
      {
        Fem=subset(b,SEX=="F")
        Male=subset(b,SEX=="M")
        if(min(c(nrow(Fem),nrow(Male)))>=N)
        {
          #FL to TL
          RLM_FL.to.TL_fem=rlm(TL~FL,data=Fem)
          R2_FL.to.TL_fem=r2(RLM_FL.to.TL_fem)
          
          RLM_FL.to.TL_male=rlm(TL~FL,data=Male)
          R2_FL.to.TL_male=r2(RLM_FL.to.TL_male)
          
          return(list(Model_FL.to.TL_fem=RLM_FL.to.TL_fem,R2_FL.to.TL_fem=R2_FL.to.TL_fem,
                      Model_FL.to.T_maleL=RLM_FL.to.TL_male,R2_FL.to.TL_male=R2_FL.to.TL_male))
        }
      }else{
        #FL to TL
        RLM_FL.to.TL=rlm(TL~FL,data=b)
        R2_FL.to.TL=r2(RLM_FL.to.TL)
        
        #       #TL to FL
        #       RLM_TL.to.FL<- rlm(FL ~ TL, data = b)
        #       R2_TL.to.FL=r2(RLM_TL.to.FL)
        
        return(list(Model_FL.to.TL=RLM_FL.to.TL,R2_FL.to.TL=R2_FL.to.TL,n=nrow(b)))
        
        #       return(list(Model_FL.to.TL=RLM_FL.to.TL,R2_FL.to.TL=R2_FL.to.TL,
        #                   Model_TL.to.FL=RLM_TL.to.FL,R2_TL.to.FL=R2_TL.to.FL))
        
      }
    }else{
      return(list(NULL))
    }
  }
}

