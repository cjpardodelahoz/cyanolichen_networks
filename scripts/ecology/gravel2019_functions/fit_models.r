# Compute the co-occurrence and the interaction models
fit_models = function(data, Enames, funC, funL, selection) {
	
	if(sum(data$Xij)!=0 ) {

		# MODEL CO-OCCURRENCE
		if (identical(funC, C2) | identical(funC, C3)) {
		   modelC = funC(data, Enames, selection)
		} else {
		   modelC = funC(data, selection)
		}

		# MODEL INTERACTIONS
		if (identical(funL, L2)) {
		   modelL = funL(data, Enames, selection)
		} else {
		   modelL = funL(data, selection)
		}
	}

	else {
		modelC = NULL
	 	modelL = NULL
	}

	return(list(i=as.character(data$IDi[1]),j=as.character(data$IDj[1]), modelC = modelC,modelL = modelL))
}