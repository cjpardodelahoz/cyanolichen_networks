# Returns the probability of observing the interaction ij 
# independently of co-occurrence and of the environment
# Havens type of interaction
# We have to trick the glm function so that it returns a probability 
# of 1 if the species do interact at least once, and 0 if they never
# interact
L0 <- function(data, selection) {
   with(as.list(data),{

        if(sum(Lij)>0) Lij_trick = numeric(length(Lij)) + 1
          else Lij_trick = numeric(length(Lij))
        model <- glm(Lij_trick ~ 1, family = "binomial")

        if(selection) {return(step(model,trace=0))}
        else {return(model)}
   })
}

# Returns the probability of observing the interaction ij when i 
# and j are present independently of the environment
L1 <- function(data, selection) {
   with(as.list(data),{
        model <- glm(Lij[Xij==1] ~ 1, family = "binomial")
        if(selection) {return(step(model,trace=0))}
        else {return(model)}
   })
}

# Returns the probability of observing the interaction ij when i 
# and j are present for each level of the environment
L2 <- function(data, Enames, selection) {
  with(as.list(data), {
    # Subset data based on the condition Xij == 1
    subLij <- Lij[Xij == 1]
    subE <- E[Xij == 1, ]
    
    # Remove categorical variables with only one level
    valid_vars <- Enames[!sapply(Enames, function(var) {
      if (is.factor(subE[[var]]) || is.character(subE[[var]])) {
        # Check if categorical and has only one level
        length(unique(subE[[var]])) <= 1
      } else {
        # Not a categorical variable, so keep it
        FALSE
      }
    })]
    
    # Update Enames and subE to reflect valid variables only
    Enames <- valid_vars
    subE <- subE[, Enames, drop = FALSE]
    
    # Create formula for the GLM
    fmla <- as.formula(paste("subLij ~", paste(Enames, collapse = "+")))
    
    # Fit the GLM and apply stepwise selection if specified
    model <- glm(fmla, family = "binomial", data = subE)
    if (selection) {
      return(step(model, trace = 0, direction = "forward"))
    } else {
      return(model)
    }
  })
}

