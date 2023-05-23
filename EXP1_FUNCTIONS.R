
################################################################   
####################### USEFUL FUNCTIONS ####################### 
################################################################ 
############## THANKS TO YOUR GOOD PAL ADRIANO :) ############## 
################################################################  


library(report)



## TEST MODEL RESIDUALS FOR NORMALITY ##

test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}


## REPORT LMER/GLMER MODEL OUTPUT ## 


output_lmer <- function(model,show_report = F) {
  print(paste("############### MODEL :",deparse(substitute(m1)),"###############",sep = " "))
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(Anova(model))
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  print("------------REPORT------------")
  if (show_report) {
    print(report(model))
  }else{ cat("Set show_report to TRUE to obtain a human readable model report")}
  #tab_model(model)
}




## SIMPLIFY MODEL IF INTERACTION NON SIGNIFICANT ##


simplify_model <- function(model) {
  # Extract the model formula
  model_formula <- formula(model)
  # Check the significance of the interaction using anova()
  anova_m1 <- as.data.frame(Anova(model))
  print(anova_m1)
  # Find the interaction term in the model formula
  # define a regular expression pattern to match the desired substring
  pattern <- "\\b\\w+\\s*\\*\\s*\\w+\\b"
  # use the sub() function to extract the first match of the pattern
  interaction_term <- sub(paste0(".*(", pattern, ").*"), "\\1", as.character(model_formula)[3])
  interaction_vars <- unlist(strsplit(interaction_term, " * "))
  anova_term <- gsub("\\s*\\*\\s*", ":", interaction_term)
  # If the Anova of the interaction is not significant, simplify the model by removing the interaction
  if (anova_m1[which(rownames(anova_m1)== anova_term),"Pr(>Chisq)"] > 0.05) {
    cat("\n#\nModel interaction NOT significant, simplify\n#\n")
    model_no_interaction_formula <-  as.formula(gsub("\\*", "+", deparse(model_formula)))
    model_no_interaction <- update(model, formula = model_no_interaction_formula)
    #print(summary(m1_no_interaction))
    print(Anova(model_no_interaction))
    return(model_no_interaction)
  } else {
    cat("\n#\nModel interaction significant, don't simplify\n#\n")
    return(model)
  }
}





## CHECK IF MODEL HAS INTERACTION ##


has_interaction <- function(model) {
  formula_str <- as.character(formula(model))
  return(any(grepl(":", formula_str) | grepl("\\*", formula_str)))
}





## COMPUTE AND STORE POST HOCS ##

# function to perform posthocs
posthoc_list <- list()
interactions_to_explore <- list()
compute_posthocs <- function(model) {
  #warning("this function has only been tested with lmer()")
  warning("How to use: \nA.Run on models without interactions. If interaction is present, run 'simplify_model' first. 
          \nB.If has_interaction= T, paste your variables naming the new var in this format 'VAR1_VAR2'
          \nC. assign the output of this function to 'posthoc_list'
          \nD.Optional: provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] to later recall the posthoc.\n\n")
  print(paste("model predictor:", paste0(row.names(Anova(model))), sep = " "))
  # check that there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) == 0) {
    print("there are no significant vars.")
  } else {
    for (SIG.VAR in row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) {
      if (grepl(":", SIG.VAR)) {
        warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function. Re-run model using pasted variables"))
        #interactions_to_explore <- c(interactions_to_explore, list(paste(GENE,GROUP,SIG.VAR,deparse(substitute(model)), sep = "-") ))
      } else {
        # check if the variable is not numeric . to do so, we need to access the dataframe from the model
        if (!is.numeric(get(gsub("\\[.*", "", as.character(model@call)[3]))[, SIG.VAR])) {
          print(paste0("Performing posthocs for the significant var: ", SIG.VAR))
          arg <- list("Tukey")
          names(arg) <- SIG.VAR
          # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
          cmp <- do.call(mcp, arg)
          posthoc_SIG.VAR <- summary(glht(model, linfct = cmp), test = adjusted("BH"))
          # Set up a compact letter display of all pair-wise comparisons
          model_means_cld <- cld(posthoc_SIG.VAR)
          # create dataframe usable with ggplot geom_text
          model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
          # add column name
          model_means_cld$newcol <- NA
          colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
          colnames(model_means_cld)[which(names(model_means_cld) == "V1")] <- "letters"
          model_means_cld[, SIG.VAR] <- row.names(model_means_cld)
          rownames(model_means_cld) <- NULL
          # if interaction term, split columns
          if (grepl("_", SIG.VAR)) {
            # Split the column into two new columns using "_"
            #SIG.VAR.STRIP <- names(model_means_cld[grep( "_",model_means_cld)])
            model_means_cld[, strsplit(SIG.VAR, "_")[[1]]] <- t(sapply(model_means_cld[,SIG.VAR], function(x) strsplit(x, "_")[[1]]))
          }
          # add to list
          posthoc_list <- c(posthoc_list, list(model_means_cld))
          if (exists("ID_model")) {
            names(posthoc_list)[length(posthoc_list)] <- paste(ID_model,SIG.VAR,  deparse(substitute(model)),sep = "-")
          } else {
            warning("if you provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] before running the function, it will be easier to later call the posthoc output for plotting")
            names(posthoc_list)[length(posthoc_list)] <- paste(SIG.VAR,deparse(substitute(model)), sep = "-")
          }
          print(paste(deparse(substitute(model)), SIG.VAR, sep = "_"))
          print(model_means_cld)
        } # SIG.VAR LOOP
      } # check if is an interaction
    } # check if numeric
    warning("call 'posthoc_list' to get posthocs")
  } # if significant vars exist
  return(posthoc_list)
}




## ADD STARS TO SIGNIFICANCE ##

add_star <- function(p) {
  if (p<0.001) {
    return('***')
  } else if (p<0.01) {
    return('**')
  } else if (p<0.05) {
    return('*')
  } else {
    return('ns')
  }
}









