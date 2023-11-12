library(rjags)
library(coda)


jags_func = function(y,n,vacc.name,placebo){
	cat("=============================\n")
	cat(paste("Vaccine:",vacc.name,"// placebo =",placebo,"\n"))
    data_obs= c(rep(1,y),rep(0,n-y))

    data   = NULL
    data$X = data_obs
    data$n = n
	
    model = "model.bug"

    jm = jags.model(model,data)

    #Burn-in
    update(jm,1000)
    chain = coda.samples(jm,c("p","y"),n.iter=10000)
    print(summary(chain))

    #Control plots
    plot(chain, col="navy") 
    mtext(paste("Vaccine:",vacc.name,"// placebo =",placebo),line=3)

    # Let"s format our chain 
    chain.df = as.data.frame(as.mcmc(chain))
  
    hist(chain.df$p, breaks = 50, prob = TRUE, col = "darkolivegreen2", xlab = "p", ylab = "f(p)", main = "Inference on p")
    mtext(paste("Vaccine:",vacc.name,"// placebo =",placebo),line=3)
    
    
    # Calcola l'intervallo di credibilità del 95%
    #Coda
    cred_interval = HPDinterval(as.mcmc(chain), prob = 0.95)

    print(cred_interval)

    

    	
}

cat("=============================\n")
cat("Jcovden // Placebo trial\n")

# Jcovden

#v stands for vaccinated; p for placebo
y_v = 116
n_v = 19630

y_p = 348
n_p = 19691


jags_func(y_p,n_p,"Jcovden",placebo="TRUE")


cat("=============================\n")
cat("Jcovden // Vaccine trial\n")

jags_func(y_v,n_v,"Jcovden",placebo="FALSE")



# Moderna/Spikevax

#v stands for vaccinated; p for placebo
y_v = 11
n_v = 14134

y_p = 185
n_p = 14073

cat("=============================\n")
cat("Moderna/Spikevax // Placebo trial\n")

jags_func(y_p,n_p,"Moderna/Spikevax",placebo="TRUE")


cat("=============================\n")
cat("Moderna/Spikevax // Vaccine trial\n")

jags_func(y_v,n_v,"Moderna/Spikevax",placebo="FALSE")


# Astrazeneca webpage no longer exists
#but I found this page: https://www.ema.europa.eu/en/news/ema-recommends-covid-19-vaccine-astrazeneca-authorisation-eu


#v stands for vaccinated; p for placebo
y_v = 64
n_v = 5258

y_p = 154
n_p = 5210

cat("=============================\n")
cat("Astrazeneca // Placebo trial\n")

jags_func(y_p,n_p,"Astrazeneca",placebo="TRUE")


cat("=============================\n")
cat("Astrazeneca // Vaccine trial\n")

jags_func(y_v,n_v,"Astrazeneca",placebo="FALSE")


