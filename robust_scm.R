# ============================================================================

# Programme: synthetic_control_monte_carlo
# Authors: Cole Dreier, Ali Uppal, Steven Yee
# Date: 07/02/2021
# Purpose: ECON 220E Monte Carlo assignment

# This programme has the following sections:
# (0)  Preliminaries
# (1)  Synthetic control function
# (2)  Randomisation inference function
# (3)  Generate data function
# (4)  Generate treatment functions
# (5)  Generate treatment application function
# (6)  Generate monte carlo simulation function
# (7)  Generate synth control chart function
# (8)  Generate size control checking function
# (9)  Generate power curve creation function
# (10) Conduct monte carlo analysis

# ============================================================================

# ===================
# (0) Preliminaries
# ===================

rm(list = ls())

# install.packages("pracma")
# install.packages("matrixcalc")
# install.packages("limSolve")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("directlabels")
library("pracma")
library("matrixcalc")
library("limSolve")
library("tidyverse")
library("ggplot2")
library("directlabels")

# ==============================
# (1) Synthetic control function
# ==============================
# This function delivers the optimal weights for the synthetic control
sc <- function(treated, control){
    # Number of units in donor pool is J0 (i.e., control units)
    J0 = dim(control)[2] 
    
    # Row vector of 1s (number of columns = units in donor pool)
    e = matrix(1, 1, J0)
    
    # Constraints that weights must sum to 1 (ex=f)
    f = 1
    
    # Square matrix (JxJ) with with ones along diagonal
    g = diag(x=1, J0, J0)
    
    # Weights are non-negative (gx>=h) (removing this improves fit, but by extrapolating!)
    h = matrix(0, J0, 1)
    
    # Least Squares with Equalities and Inequalities to solve weights
    weights = lsei(A=control, B=treated,E=e,F=f,G=g,H=h, type=2)$X
    return(list(weights=weights))
}

# ====================================
# (2) Randomisation inference function
# ====================================
# This function generates a p-value using randomisation inference
ri <- function(observed,T0,T1){
    # Number of total units (treated + control)
    J1 = dim(observed)[2]
    
    # Generate empty matrix to store the actual treated and synthetic control
    actual_synth_control = rep(NA,T1)
    actual_treated = rep(NA,T1)
    
    # Test statistics under consideration (generate empty matrices)
    RMSPE = rep(NA,J1)
    T_stat = rep(NA,J1)
    Post_treatment_avg = rep(NA,J1)
    
    # RI requires looping over all units (assigning each one the treatment status)
    for (j in 1:(J1)){
        # Making jth unit the "treated"
        treated = observed[,j] 
        
        # Excludes the jth column and uses everything else as control
        control = observed[,-j] 
        
        # Caluculate the synthetic control and weights using sc function
        weights = sc(treated[1:T0],control[1:T0,])$weights
        synth_control = control%*%weights
        
        # Calculate treatment effect (treatment takes place at T0+1)
        effect = treated - synth_control
        effect_pre = effect[1:T0]
        effect_post = effect[(T0+1):T1]
        
        # Results from each test statitic for a given unit j
        ## RMSPE test
        RMSPE[j] = sqrt(mean(effect_post^2))/sqrt(mean(effect_pre^2))
        ## T-statistic
        average_treat = mean(effect[(T0+1):T1])
        stdev_treat = std(effect[(T0+1):T1])
        T_stat[j] = abs((average_treat/(T1-T0))/(stdev_treat/sqrt(T1-T0)))
        ## Post-treatment average
        Post_treatment_avg[j] = mean(abs(effect[(T0+1):T1]))
        if(j==1){
            actual_synth_control = synth_control
            actual_treated = treated
        }
    }
    
    # Pvalue for each test stat is proportion of placebos with a greater test stat
    # than actual treated unit which is the first entry in the matrix
    pvalue_RMPSE = mean(RMSPE>=RMSPE[1])
    pvalue_tstat = mean(T_stat >= T_stat[1])
    pvalue_post = mean(Post_treatment_avg >= Post_treatment_avg[1])
    
    return(list(pvalue_RMSPE=pvalue_RMPSE, pvalue_tstat=pvalue_tstat, pvalue_post=pvalue_post, actual_treated=actual_treated, actual_synth_control=actual_synth_control))
}

# ==========================
# (3) Generate data function
# ==========================
# This function uses different DGPs to generate the data.
gen_data<- function(DGP, J1, T1){
    # Generate an empty data matrix where rows = time periods, cols = units
    data_mat = matrix(NA,T1,J1)
    
    # if DGP=5, draw from normal dist with twfe. 
    if (DGP == 5){
        # Creating T1xJ1 iid shocks
        eps = matrix(rnorm(J1*T1,0,5), T1, J1)
        #eps = matrix(rexp(J1*T1),T1,J1)
        
        # Defining time fixed effects
        time_fe = matrix(rnorm(T1,0,1), T1, 1)
        
        # Defining unit fixed effects
        unit_fe = matrix(rnorm(T1,0,1), 1, J1)
        
        # Generate common factors
        time_fe_mat = repmat(time_fe, 1, J1)
        unit_fe = repmat(unit_fe, T1, 1)
        
        # Generate untreated outcomes
        data_mat = unit_fe + time_fe_mat + eps
    }
    
    
    # if DGP=4, draw from exponention dist. 
    if (DGP == 4){
        data_mat = matrix(rexp(J1*T1),T1,J1)
    }
    
    # if DGP=3, draw from normal dist with twfe and unit_fe of treatment unit < unit_fe of all others. 
    else if(DGP == 3){
        # Creating T1xJ1 iid shocks
        eps = matrix(rnorm(J1*T1,0,1), T1, J1)
        
        # Defining time fixed effects
        time_fe = matrix(rnorm(T1,0,1), T1, 1)
        
        # Defining unit fixed effects
        unit_fe = matrix(NA, 1, J1)
        for (i in 2:J1){
            unit_fe[,i] = rnorm(1,0,1)
        }
        # Set treatment unit's fe as 10% lower than the minimum of control unit fe
        unit_fe[,1] = -1.1*abs(max(unit_fe[2:J1])) 
        
        # Generate common factors
        time_fe_mat = repmat(time_fe, 1, J1)
        unit_fe = repmat(unit_fe, T1, 1)
        
        # Generate untreated outcomes
        data_mat = unit_fe + time_fe_mat + eps
    }
    
    # if DGP=2, draw from normal dist with twfe and unit_fe of treatment unit is the average of unit_fe of all others. 
    else if(DGP == 2){
        # Creating T1xJ1 iid shocks
        eps = matrix(rnorm(J1*T1,0,1), T1, J1)
        
        # Defining time fixed effects
        time_fe = matrix(rnorm(T1,0,1), T1, 1)
        
        # Defining unit fixed effects
        unit_fe = matrix(NA, 1, J1)
        for (i in 2:J1){
            unit_fe[,i] = rnorm(1,0,1)
        }
        # Set treated unit to be in the middle of the distribution of control unit fe
        unit_fe[,1] = mean(unit_fe[2:J1]) 
        
        # Generate common factors
        time_fe_mat = repmat(time_fe, 1, J1)
        unit_fe = repmat(unit_fe, T1, 1)
        
        # Generate untreated outcomes
        data_mat = unit_fe + time_fe_mat + eps
    }
    
    # if DGP=1, draw from normal dist only. 
    else if (DGP == 1){
        data_mat = matrix(rnorm(J1*T1),T1,J1)
    }
    return(data_mat)
}

# ================================
# (4) Generate treatment functions
# ================================
# This function generates a vector of a constant treatment shifter
gen_treatment <- function(DGP, lambda, T0, T1){
    # Treatment vector is zero until treatment period
    treatment_vec = rep(0, T1)
    
    # Treatment vector is lambda every period after treatment period
    treatment_vec[(T0+1):T1] = lambda
    return(treatment_vec)
}

# This function generates a vector of a time-varying treatment shifter
gen_treatment_varying <- function(DGP, lambda, T0, T1){
    # Treatment vector is zero until treatment period
    treatment_vec = rep(0, T1)
    
    # Treatment vector is lambda every period after treatment period
    treatment_vec[(T0+1):T1] = lambda
    
    # Set treatment vector to be increasing linearly over time
    time_offset = linspace(-T0+1, T1-T0,T1)
    treatment_vec = hadamard.prod(treatment_vec, time_offset)
    return(treatment_vec)
}

# ===========================================
# (5) Generate treatment application function
# ===========================================
# This function applies the treatment shifter to the treated data
apply_treatment <- function(DGP, observed, treatment_vec){
    treated_data = observed
    
    # Add treatment vector to 1st column of untreated values. Thus, 1st column is treated unit
    treated_data[,1] = treated_data[,1] + treatment_vec 
    return(treated_data)
}

# ============================================
# (6) Generate monte carlo simulation function
# ============================================
# This function runs a monte carlo simulation
simulate <- function(DGP,sims,lambda_vals,lambda_start,lambda_end,varying,T0,T1,J0,J1){
    
    # Generate all possible values of treatment shifter (lambda)
    lambda_seq = linspace(lambda_start, lambda_end, lambda_vals+1)
    
    # Generate matrices for the different test stats
    pvalue_RMSPE_mat = matrix(NA,sims,lambda_vals+1)
    pvalue_tstat_mat = matrix(NA,sims,lambda_vals+1)
    pvalue_post_mat = matrix(NA,sims,lambda_vals+1)
    
    # Generate matrices for the actual treated and actual synthetic control
    synth_sum = matrix(0, lambda_vals+1, T1)
    treated_sum = matrix(0, lambda_vals+1, T1)
    
    # Loop first over the different treatment shifters, then over the no. of simulations
    for (iter1 in 1:(lambda_vals+1)){
        lambda_test = lambda_seq[iter1]
        for (iter2 in 1:sims){
            
            # for each simulation generate a data matrix using the data function
            data_mat = gen_data(DGP, J1, T1)
            
            # Determine if treatment shifter is time-varying; choose relevant treatment function
            if (varying){
                treatment_vec = gen_treatment_varying(DGP, lambda_test, T0, T1)
            }
            else{ 
                treatment_vec = gen_treatment(DGP, lambda_test, T0, T1)
            }
            
            # Observed data is what happens to the data matrix when applying treatment function
            observed = apply_treatment(DGP, data_mat, treatment_vec)
            
            # Apply randomisation inference function to generate p-values
            results = ri(observed,T0,T1)
            pvalue_RMSPE_mat[iter2, iter1] = results$pvalue_RMSPE
            pvalue_tstat_mat[iter2, iter1] = results$pvalue_tstat
            pvalue_post_mat[iter2, iter1] = results$pvalue_post
            
            # For a given lambda, sum across all simulations of actual synth & actual treatment 
            synth_sum[iter1,] = synth_sum[iter1,] + results$actual_synth_control
            treated_sum[iter1,] = treated_sum[iter1,] + results$actual_treated
        }
    }
    
    # For a given lambda and each time period, calculate the avg actual synth & treatment  
    synth_avg = synth_sum / sims
    treated_avg = treated_sum / sims
    return(list(pvalue_RMSPE_mat=pvalue_RMSPE_mat, pvalue_tstat_mat=pvalue_tstat_mat,pvalue_post_mat=pvalue_post_mat, synth_avg=synth_avg, treated_avg=treated_avg))
}

# =========================================
# (7) Generate synth control chart function
# =========================================
# This function creates the typical synthetic control charts
synth_chart <-function(varying, treated_avg, synth_avg, T1, lambda_end, lambda_val,DGP){
    if (lambda_val==0){
        index=1
    }
    else if (lambda_val==0.5){
        index=6
    }
    else {
        index=11
    }
    # Make dataframe
    df = data.frame(treatment_group=c(treated_avg[index,]),synthetic_control=c(synth_avg[index,]),Time=linspace(1,T1,T1))
    
    # Convert dataframe from wide to long
    df2 = gather(df,Group,Outcome,treatment_group:synthetic_control,factor_key=TRUE)
    
    # Different DGPs need different axes
    if (varying){
        ymin = -2.3
        ymax = 10*lambda_end+0.2
    }
    else {
        ymin = -2.3
        ymax = lambda_end+0.2
    }
    
    # Create ggplot synth chart
    ggplot(data=df2,aes(x=Time,y=Outcome, group=Group,colour=Group))+
        geom_line() +
        scale_colour_discrete(guide = 'none') +
        scale_x_continuous(limits = c(0, 25), expand = expansion(mult = c(0, 0.1), add = c(0, 5))) +
        scale_y_continuous(limits = c(ymin, ymax)) +
        geom_dl(aes(label = Group),  method = list(dl.trans(x = x + .1), "last.bumpup",cex = 0.65))
    
    # Saving the chart
    if (varying){
        varying1 = "TRUE"
    }
    else {
        varying1 = "FALSE"
    }
    filename=paste("synthcontrol_","DGP_",toString(DGP),"_",varying1 ,"_",toString(lambda_val),".png",sep="")
    ggsave(filename=filename,plot=last_plot())
}

# ==================================
# (8) Size control checking function
# ==================================
# This function creates a chart looking at size control across varying test sizes
size_control <-function(DGP, varying,lambda_vals,pvalue_RMSPE_mat, pvalue_tstat_mat, pvalue_post_mat,size_vals){
    size_seq = linspace(0,1,size_vals+1)
    pvalue_plot = matrix(NA,size_vals+1,4)
    
    for (i in 0:size_vals+1){
        pvalue_plot[i,1]=size_seq[i]
        pvalue_plot[i,2]=mean(pvalue_RMSPE_mat[,1] <= size_seq[i])
        pvalue_plot[i,3]=mean(pvalue_tstat_mat[,1] <= size_seq[i])
        pvalue_plot[i,4]=mean(pvalue_post_mat[,1] <= size_seq[i])
    }
    
    # Make dataframe
    df = data.frame(Size=c(pvalue_plot[,1]),RMSPE=c(pvalue_plot[,2]), Tstat=c(pvalue_plot[,3]), Post=c(pvalue_plot[,4]))
    
    # Convert dataframe from wide to long
    df2 = gather(df,Test,Result,RMSPE:Post,factor_key=TRUE)
    
    # Create ggplot synth chart
    ggplot(data=df2,aes(x=Size,y=Result, group=Test,colour=Test))+
        geom_line() +
        scale_colour_discrete(guide = 'none') +
        scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0), add = c(0, 0.24))) +
        scale_y_continuous(limits = c(0, 1.26)) +
        geom_dl(aes(label = Test),  method = list(dl.trans(x = x + .1), "last.bumpup",cex = 0.8))
    
    # Saving the chart
    filename=paste("sizecontrol_","DGP_",toString(DGP),".png",sep="")
    ggsave(filename=filename,plot=last_plot())
}

# ==========================================
# (9) Generate power curve creation function
# ==========================================
# This function creates a chart looking at power for different effect sizes
gen_power_curve <- function(DGP, varying,lambda_vals,lambda_start,lambda_end, pvalue_RMSPE_mat, pvalue_tstat_mat, pvalue_post_mat, size){
    lambda_seq = linspace(lambda_start, lambda_end, lambda_vals+1)
    pvalue_plot = matrix(NA,lambda_vals+1,4)
    
    for (i in 0:lambda_vals+1){
        pvalue_plot[i,1]=lambda_seq[i]
        pvalue_plot[i,2]=mean(pvalue_RMSPE_mat[,i] <= size)
        pvalue_plot[i,3]=mean(pvalue_tstat_mat[,i] <= size)
        pvalue_plot[i,4]=mean(pvalue_post_mat[,i] <= size)
    }
    
    # Make dataframe
    df = data.frame(Effect=c(pvalue_plot[,1]),RMSPE=c(pvalue_plot[,2]), Tstat=c(pvalue_plot[,3]), Post=c(pvalue_plot[,4]))
    
    # Convert dataframe from wide to long
    df2 = gather(df,Test,Power,RMSPE:Post,factor_key=TRUE)
    
    # Create ggplot synth chart
    ggplot(data=df2,aes(x=Effect,y=Power, group=Test,colour=Test))+
        geom_line() +
        scale_colour_discrete(guide = 'none') +
        scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0), add = c(0, 0.24))) +
        scale_y_continuous(limits = c(0, 1.26)) +
        geom_dl(aes(label = Test),  method = list(dl.trans(x = x + .1), "last.bumpup",cex = 0.8))
    
    # Saving the chart
    if (varying){
        varying1 = "TRUE"
    }
    else {
        varying1 = "FALSE"
    }
    filename=paste("powercurve_DGP_",toString(DGP),"_",varying1,".png",sep="")
    ggsave(filename=filename,plot=last_plot())
}

# =================================
# (10) Conduct monte carlo analysis
# =================================
# Set simulation parameters
T1  = 25
T0 = 15
J0 = 19
J1 = J0 + 1
sims = 1000
lambda_vals = 10
lambda_start = 0
lambda_end = 1
size_vals = 10
size = 0.1

# Run simulation
varying = FALSE
for (j in 1:2){
    # If j=1, lambda is constant. Now test each DGP (1,2,3)
    for (i in 1:3){
        simulation = simulate(DGP=i,sims,lambda_vals,lambda_start,lambda_end, varying, T0,T1,J0,J1)
        
        # Create synth charts across each value of lambda (note here this only works for default lambda values)
        for (k in c(0,0.5,1)){
            synth_chart(varying,simulation$treated_avg,simulation$synth_avg,T1,lambda_end,lambda_val=k,DGP=i)
        }
        
        # Check size control only if the initial lambda value = 0 (only do for j=1 as it is same for j=2)
        if (j==1 && linspace(lambda_start, lambda_end, lambda_vals+1)[1]==0){
            size_control(DGP=i,varying,lambda_vals,simulation$pvalue_RMSPE_mat, simulation$pvalue_tstat_mat, simulation$pvalue_post_mat,size_vals)
        }
        
        # Create power curve for each DGP
        gen_power_curve(DGP=i, varying, lambda_vals,lambda_start,lambda_end, simulation$pvalue_RMSPE_mat, simulation$pvalue_tstat_mat, simulation$pvalue_post_mat, size)
    }
    
    # Repeat above with lambda varying
    varying=TRUE
}