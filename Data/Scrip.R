########################################################################
#######################################################################
####### Selection Index ###############################################

######################################################################
### Load library
library(tidyverse)
library(lme4)

### Load file
df <- read.csv("C1_PSI_05_Pheno.csv", header = TRUE)

### calculate the  mean by each entry (genotype)
entry_mean<- df %>%
  group_by(entry) %>% 
    summarise_at(vars(starts_with("T")),~mean(.)) %>% 
      select(-T4)

#### names of each traits that will calculate  the G y P matrix
Traits <- c("T1","T2","T3")
#### functions to obtain the valus 
co_gp<- function(Traits,Entry,Rep,data){
  
  traits <- Traits
  geno<-as.factor(Entry)
  rep <-as.factor(Rep)
  ###########################################  
  nrep <- length(levels(rep))
  leng_traits<- length(traits) # number of traits
  
  ###########################################
  
  # Creating the matrix  G and P
  
  
  G <- matrix(nrow = leng_traits, ncol = leng_traits) 
  P <- matrix(nrow = leng_traits, ncol = leng_traits) 
  # Estimation of  variance each traits.
  
  for (i in 1: leng_traits) {
    y <- data[,traits[i]]
    fm <- lmer(y ~ (1|geno)+(1|rep))
    vc <- VarCorr(fm,comp="Variance")
    G[i, i] <- vc$geno[1]
    P[i, i] <- vc$geno[1] + attr(vc, "sc")^2/nrep 
  }
  
  ####Sum de each variable  and  estimation of  variance components
  ###Example X1+X2=X1X2, X1+X3=X1X3 etc
  
  for (i in 1:( leng_traits - 1)) {
    for (j in (i + 1): leng_traits) {
      y<-data[,traits[i]] + data[,traits[j]]
      fm <-lmer(y ~ (1|  geno) + (1| rep))
      varcor <-VarCorr(fm) 
      G[i, j] <- G[j, i] <- (varcor$geno[1] - G[i, i] - G[j, j]) / 2
      P[i, j] <- P[j, i] <- (varcor$geno[1] + attr(varcor, "sc")^2 / nrep - P[i, i] - P[j, j]) / 2
      
    }
  }
  
  ####################Estimation of correlation and covariance
  diag_G <- diag(diag(G)^{-0.5},  leng_traits,  leng_traits)
  diag_P<- diag(diag(P)^{-0.5},  leng_traits,  leng_traits)
  GC <- diag_G %*% G %*% diag_G# Genotypic correlation matrix
  PC <- diag_P %*% P %*% diag_P # Phenotypic correlation matrix
  ####################Names of  matrix    
  row.names(G) <- Traits
  colnames(G) <- Traits
  row.names(P) <- Traits
  colnames(P) <- Traits
  row.names(GC) <- Traits
  colnames(GC) <- Traits
  row.names(PC) <- Traits
  colnames(PC) <- Traits
  
  G <- round(G,2)
  P <- round(P,2)
  GC <- round(GC,2)
  PC <- round(PC,2)
  # results
  results<- list(Genetic_Cov = G, Pheno_Cov = P, Genetic_Cor = GC, Pheno_Cor = PC)
}
#### run the function
df_matrix<- co_gp(Traits,df$entry,df$rep,df)
#### obtain the matrix P and G 
G<- vc$Genetic_Cov 
P <- vc$Pheno_Cov
#### Inversa of the phenotypic matrix
P1 <- solve(P)
#### vector of the weight
w<- matrix(c(1,-1,1), nrow = 3,ncol = 1)
#### obatin the b valus b = w`GP
b<- t(w)%*%(G%*%P1)
#### caculate the value I for each genotype
I<- entry_mean %>% 
  mutate(t1=T1*b[1],
         t2=T2*b[2],
         t3=T3*b[3],
         I=t1+t2+t3) %>% 
    select(-t1,-t2,-t3)
I
#### results of the table 2.3
table2.3 <- I%>% 
              filter(entry<21) %>% 
                    arrange(-I) %>%     
                      select(entry,T2,T2,T3,I)

###### selecting the bests genotype
table2.4 <- I%>% 
              arrange(-I) 



table2.4 %>% 
  arrange(-I) %>% 
      top_n(25) %>% 
  summarise(mean_t1=mean(T1),mean_t2=mean(T2),mean_t3=mean(T3),
            mean_i=mean(I))
### expected gain for 5 % of selecting
E<- 2.063*b%*%G/9.312
#### variance   and standar desviation of the estimated selecting index 
V<- b%*%(P%*%t(b))
<- t(w)%*%(G%*%w)


I%*%(G%*%t(I))
I%*%(P%*%t(I))




sd_V<- sqrt(V)
R <- 2.063*(sd_V)


9.31/10.41

###########################################################
###########################################################
############## BASE INDEX
y <- as.matrix(df[4:6])



v <- qnorm(0.5) 
k < dnorm(v)/q
 
qnorm()