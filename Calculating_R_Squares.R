###########################################################################################################################
# File Name : Calculating_R_Squares.R                                                                                     #
# Function : This script is based on the article of S.Nakagawa and H.Schielzeth published in 2013 in Methods in Ecology   #
# and Evolution : A general and simple method for obtaining R square from generalized linear mixed-effects models. And    #
# on "The coefficient of determination R2 and intra-class correlation from generalized linear mixed-effects models        #
# revisited and expanded", by the same authors and P.Johnson, published in 2017.                                          #
# This script allows to compute marginal and conditional R square for linear and generalized mixed models                 #
# For more informations : https://www.researchgate.net/publication/320387997_Worked_examples_Appendix_S6                  #
# Author : Emilie Fleurot                                                                                                 #
# Github : emilie-fleurot                                                                                                 #
# Date : 02/01/2021                                                                                                       #
###########################################################################################################################

Calculating_R_Squares<-function(model,model_nul)
{
  # Packages
  library(lme4)

  if(model@call[[1]]== "glmer")
  {
    if(model@call$family == "binomial")
    {
      # First we want to compute the marginal R-square i.e the variance part explained by fixed effects
      VarF <- var(as.vector(fixef(model) %*% t(model@pp$X))) # Estimating VarF (Variance due to Fixed effects) can,in principle, 
      # be carried out by predicting fitted values based on the fixed effects alone (equivalent to multiplying the design matrix of the fixed effects 
      # (=t(model@pp$X)) with the vector of fixed effect estimates (=as.vector(fixef(model))), followed by calculating the variance of those 
      # fitted values (var))
      
      # VarCorr() extracts variance components with a correction for binomial distribution : (pi^2)/3
      # VarCorr()$plot extract the variance of the plot effect
      # We store the variance explained by random effects in model in VarR and in the model model_nul in VarR_nul :
      VarR<-sum(as.numeric(VarCorr(model)))
      VarR_nul<-sum(as.numeric(VarCorr(model_nul)))

      Marginal_R_Square_binom <- VarF/(VarF + VarR + (pi^2)/3)
      
      # compute the conditional R-square (i.e percent of variance explained by fixed and random effects) with a correction for binomial 
      # distribution : (pi^2)/3
      Conditional_R_Square_binom <- (VarF + VarR)/(VarF + VarR + (pi^2)/3)
      
      # computing the percent of explained variance
      # for the plot slope
      1-((VarR)^2/(VarR_nul)^2)
      #for the residuals
      1-(var(residuals(model))/var(residuals(model_nul)))
      
      print(c("Marginal_R_Square_binom",Marginal_R_Square_binom))
      print(c("Conditional_R_Square_binom",Conditional_R_Square_binom))
    }
    
    if(model@call$family == "poisson")
    {
      # First we want to compute the marginal R-square i.e the variance part explained by fixed effects
      VarF <- var(as.vector(fixef(model) %*% t(model@pp$X))) # Estimating VarF (Variance due to Fixed effects) can,in principle, 
      # be carried out by predicting fitted values based on the fixed effects alone (equivalent to multiplying the design matrix of the fixed effects 
      # (=t(model@pp$X)) with the vector of fixed effect estimates (=as.vector(fixef(model))), followed by calculating the variance of those 
      # fitted values (var))
      
      # VarCorr() extracts variance components with a correction for poisson distribution : log(1 + 1/exp(as.numeric(fixef(model_nul))))
      # VarCorr()$plot extract the variance of the plot effect
      # We store the variance explained by random effects in model in VarR and in the model model_nul in VarR_nul :
      VarR<-sum(as.numeric(VarCorr(model)))
      VarR_nul<-sum(as.numeric(VarCorr(model_nul)))
      
      Marginal_R_Square_poisson <- VarF/(VarF + VarR + log(1 + 1/exp(as.numeric(fixef(model_nul)))))
      
      # compute the conditional R-square (i.e percent of variance explained by fixed and random effects) with a correction for poisson 
      # distribution : log(1 + 1/exp(as.numeric(fixef(model_nul))))
      Conditional_R_Square_poisson <- (VarF + VarR)/(VarF + VarR + log(1 + 1/exp(as.numeric(fixef(model_nul)))))
      
      # computing the percent of explained variance
      # for the plot slope
      1-((VarR)^2/(VarR_nul)^2)
      #for the residuals
      1-(var(residuals(model))/var(residuals(model_nul))) 
      
      print(c("Marginal_R_Square_poisson",Marginal_R_Square_poisson))
      print(c("Conditional_R_Square_poisson",Conditional_R_Square_poisson))
    }
    
    if(length(grep("MASS::negative.binomial",model@call$family)) == 1) # If the distribution used in the model is a negative.binomial
    {   # With the function grep() we make sure that the family name of the model contains the pattern "MASS::negative.binomial"
      # First we want to compute the marginal R-square i.e the variance part explained by fixed effects
      VarF <- var(as.vector(fixef(model) %*% t(model@pp$X))) # Estimating VarF (Variance due to Fixed effects) can,in principle, 
      # be carried out by predicting fitted values based on the fixed effects alone (equivalent to multiplying the design matrix of the fixed effects 
      # (=t(model@pp$X)) with the vector of fixed effect estimates (=as.vector(fixef(model))), followed by calculating the variance of those 
      # fitted values (var))
      
      # VarCorr() extracts variance components with a correction for negative binomial distribution : VarOtF
      # VarCorr()$plot extract the variance of the plot effect
      # We store the variance explained by random effects in model in VarR and in the model model_nul in VarR_nul :
      VarR<-sum(as.numeric(VarCorr(model)))
      VarR_nul<-sum(as.numeric(VarCorr(model_nul)))
      
      # getting the observation-level variance Null model
      thetaN <- getME(model_nul, "glmer.nb.theta")
      lambda <- as.numeric(exp(fixef(model_nul) + 0.5 * VarR))
      
      # lambda2 <- mean(DataAll$Parasite)
      VarOdN <- 1/lambda + 1/thetaN # the delta method
      VarOlN <- log(1+ (1/lambda) + (1/thetaN)) # log-normal approximation
      VarOtN <- trigamma((1/lambda + 1/thetaN)^(-1)) # trigamma function
      # comparing the three
      c(VarOdN = VarOdN, VarOlN = VarOlN, VarOtN = VarOtN)
      ## VarOdN VarOlN VarOtN
      ## 0.6520918 0.5020423 0.9077864
      # Full model
      thetaF <- getME(model, "glmer.nb.theta")
      VarOdF <- 1/lambda + 1/thetaF # the delta method
      VarOlF <- log(1+ (1/lambda) + (1/thetaF)) # log-normal approximation
      VarOtF <- trigamma((1/lambda + 1/thetaF)^(-1)) # trigamma function
      # comparing the three
      c(VarOdF = VarOdF, VarOlF = VarOlF, VarOtF = VarOtF)
      
      Marginal_R_Square_nbinom <- VarF/(VarF + VarR + VarOtF)
      
      # compute the conditional R-square (i.e percent of variance explained by fixed and random effects) with a correction for negative binomial 
      # distribution : VarOtF
      Conditional_R_Square_nbinom <- (VarF + VarR)/(VarF + VarR + VarOtF)
      
      # computing the percent of explained variance
      # for the plot slope
      1-((VarR)^2/(VarR_nul)^2)
      #for the residuals
      1-(var(residuals(model))/var(residuals(model_nul)))
      
      print(c("Marginal_R_Square_nbinom",Marginal_R_Square_nbinom))
      print(c("Conditional_R_Square_nbinom",Conditional_R_Square_nbinom))
    }
  }
  
  if(model@call[[1]]== "lmer")
  {
    # First we want to compute the marginal R-square i.e the variance part explained by fixed effects
    VarF <- var(as.vector(fixef(model) %*% t(model@pp$X))) # Estimating VarF (Variance due to Fixed effects) can,in principle, 
    # be carried out by predicting fitted values based on the fixed effects alone (equivalent to multiplying the design matrix of the fixed effects 
    # (=t(model@pp$X)) with the vector of fixed effect estimates (=as.vector(fixef(model))), followed by calculating the variance of those 
    # fitted values (var))
    
    # VarCorr() extracts variance components with a correction for gaussian distribution : attr(VarCorr(model), "sc")^2
    # attr(VarCorr(lmer.model),’sc’)^2 extracts the residual variance, VarCorr()$plot extract the variance of the plot effect
    # We store the variance explained by random effects in model in VarR and in the model model_nul in VarR_nul :
    VarR<-sum(as.numeric(VarCorr(model)))
    VarR_nul<-sum(as.numeric(VarCorr(model_nul)))
    
    Marginal_R_Square_lmer <- VarF/(VarF + VarR + (attr(VarCorr(model), "sc")^2))
    
    # compute the conditional R-square (i.e percent of variance explained by fixed and random effects) with a correction for gaussian 
    # distribution : attr(VarCorr(model), "sc")^2
    Conditional_R_Square_lmer <- (VarF + VarR)/(VarF + VarR + (attr(VarCorr(model), "sc")^2))
    
    # computing the percent of explained variance
    # for the plot slope
    1-((VarR)^2/(VarR_nul)^2)
    #for the residuals
    1-(var(residuals(model))/var(residuals(model_nul)))
    
    print(c("Marginal_R_Square_lmer",Marginal_R_Square_lmer))
    print(c("Conditional_R_Square_lmer",Conditional_R_Square_lmer))
  }
  
  else
  {
  }
}

