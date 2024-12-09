fit_all_pairs <- function(DF_split, Enames) {
    results <- lapply(seq_along(DF_split), function(pair_index) {
        data <- DF_split[[pair_index]]
        
        # Pick the model
        models_C2_L0 <- fit_models(data, selection = FALSE, funC = C2, funL = L0, Enames)
        models_C2_L1 <- fit_models(data, selection = FALSE, funC = C2, funL = L1, Enames)
        models_C2_L2 <- fit_models(data, selection = FALSE, funC = C2, funL = L2, Enames)
        models_C0_L2 <- fit_models(data, selection = FALSE, funC = C0, funL = L2, Enames)
        models_C1_L2 <- fit_models(data, selection = FALSE, funC = C1, funL = L2, Enames)
        models_C3_L2 <- fit_models(data, selection = FALSE, funC = C3, funL = L2, Enames)
        
        # Compute the LL
        LL_C2_L0 <- get_LL(models_C2_L0, data)
        LL_C2_L1 <- get_LL(models_C2_L1, data)
        LL_C2_L2 <- get_LL(models_C2_L2, data)
        LL_C0_L2 <- get_LL(models_C0_L2, data)
        LL_C1_L2 <- get_LL(models_C1_L2, data)
        LL_C3_L2 <- get_LL(models_C3_L2, data)
        
        # Collect the results
        LL <- c(
            LL_C2_L0[1],
            LL_C2_L1[1],
            LL_C2_L2[1],
            LL_C0_L2[1],
            LL_C1_L2[1],
            LL_C3_L2[1]	
        )
        
        npars <- c(
            LL_C2_L0[2],
            LL_C2_L1[2],
            LL_C2_L2[2],
            LL_C0_L2[2],
            LL_C1_L2[2],
            LL_C3_L2[2]	
        )
        
        AIC <- -2 * LL + 2 * npars
        
        # Put results in a data frame
        result <- data.frame(
            LL = LL,
            npars = npars,
            AIC = AIC,
            model = c("C2_L0", "C2_L1", "C2_L2", "C0_L2", "C1_L2", "C3_L2"),
            mycobiont = sapply(strsplit(names(DF_split[pair_index]), " "), "[", 1),
            nostoc = sapply(strsplit(names(DF_split[pair_index]), " "), "[", 2)
        )

    })
    
    # Combine results for all pairs
    do.call(rbind, results)
}