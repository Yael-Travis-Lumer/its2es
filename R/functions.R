#' ITS analysis for count outcomes with no seasonality adjustment
#'
#' \code{its_poisson_wo_seas} fits a Poisson regression model to an ITS, and returns the model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param offset_name either a string indicating the name of the offset column in the data, or NULL. Default value is NULL
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param over_dispersion Logical - indicating whether a quasi-Poisson model should be used to account for over dispersion (TRUE), or not (FLASE). Default value is FALSE.
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted Poisson regression model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_poisson_wo_seas(data=data,form=form,offset_name="labour", time_name = "time",intervention_start_ind=intervention_start_ind,over_dispersion=TRUE, impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_poisson_wo_seas <- function(data,form, offset_name=NULL,time_name, intervention_start_ind,over_dispersion=FALSE, impact_model="full", counterfactual="FALSE"){
  pred <- predC <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if(!(is_formula(form))){
    stop("Please make sure that form is a formula object")
  }
  if (!(time_name %in% colnames(data))){
    stop("Please make sure that time_name belongs to colnames(data)")
  }
  if(!all(data[[time_name]]==(1:nrow(data)))){
    data[[time_name]]=1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(is.null(offset_name))){
    if (!(offset_name %in% colnames(data))){
      stop("Please make sure that offset_name belongs to colnames(data)")
    }
    if (offset_name %in% covariates | paste0("offset(log(",offset_name,"))") %in% covariates | paste0("offset(",offset_name,")") %in% covariates | paste0("log(",offset_name,")") %in% covariates){
      stop("The offset term should not be included in the formula. Please supply a formula without the offset, and add the name of the offset column using the argument offset_name.")
    }
  }
  if(!is.logical(over_dispersion)){
    stop("over_dispersion must be either TRUE or FALSE")
  }
  if (!(impact_model %in% c("full","level","slope"))){
    stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  }
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  }
  if (!all(covariates %in% colnames(data))){
    stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  data$indicator <- 0
  data$indicator[intervention_start_ind:n] <- 1
  data$shifted_time <- data[[time_name]]-intervention_start_ind
  # update formula
  if (impact_model=="full"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator + indicator:shifted_time")))
    }
  } else if (impact_model=="level"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator")))
    }
  } else if (impact_model=="slope"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name,"+ indicator:shifted_time")))
    }
  }
  if (!is.null(offset_name)){
    offset_name <- sym(offset_name)
    form_update <-update(form_update, as.formula(paste0("~ . + offset(log(",offset_name,"))")))
  }
  # fit model
  if(!isTRUE(over_dispersion)){
    set.seed(1)
    model <- glm(form_update,data,family=poisson)
  } else {
    set.seed(1)
    model <- glm(form_update,data,family=quasipoisson)
  }

  t_star <- which(data$indicator==1, arr.ind = TRUE)[1] #intervention_start_ind
  b <- (n-t_star)/2

  cov_mat <- vcov(model)
  if (impact_model=="full"){
    beta_intervention <- model$coefficients["indicator"]
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    RR <- exp(beta_intervention+beta_interaction*b)
    var_lin <- cov_mat["indicator","indicator"]+b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]+b*2*cov_mat["indicator","indicator:shifted_time"]
    lin_est <- beta_intervention+b*beta_interaction
  } else if (impact_model=="level"){
    beta_intervention <- model$coefficients["indicator"]
    RR <- exp(beta_intervention)
    var_lin <- cov_mat["indicator","indicator"]
    lin_est <- beta_intervention
  }
  else if (impact_model=="slope"){
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    RR <- exp(beta_interaction*b)
    var_lin <- b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]#+b*2*cov_mat["indicator","indicator:shifted_time"]
    lin_est <- b*beta_interaction
  }

  CI_lin <- c(lin_est-1.96*sqrt(var_lin),lin_est+1.96*sqrt(var_lin))
  CI_RR <- exp(CI_lin)
  TS <- lin_est/sqrt(var_lin)
  p_value <- round(2*pnorm(abs(TS),lower.tail = FALSE),2)
  ret_vec <-c(RR,CI_RR,p_value)
  names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI","P-value")
  s <- summary(model)
  print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec #or model[["RR"]] <- ret_vec and then return(model)

  #add predictions
  predA <- predict(model,type="response")
  data$pred <- predA#*100/data[[offset_name]]
  if (isTRUE(counterfactual)){
    # Compute counterfactual values
    new_data <- data
    new_data$indicator <- 0
    pred_counter <- predict(model,type="response",newdata = new_data)
    data$predC <- pred_counter#*100/data[[offset_name]]
    data$predC[data$indicator!=1] <- NA
  }
  data <- as_tibble(data)
  # #draw
  # Group <- levels( factor(c("Counterfactual","Fitted values")))
  # response <- sym(response)
  # p2 <- ggplot(data = data , aes(y = !!response*100/data[[offset_name]], x = dt)) +
  #   geom_point()  + ylab(y_lab) + xlab("Date") +
  #   geom_line(aes(y = predC, x = dt, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth)+
  #   geom_line(aes(y = pred, x = dt, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
  #   scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
  #   scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
  #   theme(legend.key=element_blank())
  # p2
  return(list(model=model, model_summary=s, data=data))
}

#' ITS analysis for count outcomes with seasonal adjustment via Fourier terms
#' \code{its_poisson_fourier} fits a Poisson regression model adjusted to seasonality, and returns the model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param offset_name either a string indicating the name of the offset column in the data, or NULL. Default value is NULL
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param over_dispersion Logical - stating whether an over-dispersed Poisson model should be used (TRUE), or not (FALSE). Default is FALSE and then a regular Poisson model is used. If TRUE then a quasi-Poisson model is used instead.
#' @param freq A positive integer describing the frequency of the time series.
#' @param keep_significant_fourier Logical - indicating whether only the significant Fourier terms should be considered. Default is TRUE and then the model is fitted twice; once to obtain the significant Fourier terms, and second time keeping only the significant Fourier terms. If FALSE, then all the Fourier terms are used.
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted Poisson regression model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_poisson_fourier(data=data,form=form,offset_name="labour", time_name = "time",intervention_start_ind=intervention_start_ind,over_dispersion=TRUE, freq=12, keep_significant_fourier=TRUE, impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_poisson_fourier <- function(data, form, offset_name=NULL,time_name, intervention_start_ind,over_dispersion=FALSE,freq, keep_significant_fourier=TRUE,impact_model="full", counterfactual="FALSE"){
  pred <- predC <- NULL
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if(!(is_formula(form))){
    stop("Please make sure that form is a formula object")
  }
  if (!(time_name %in% colnames(data))){
    stop("Please make sure that time_name belongs to colnames(data)")
  }
  if(!all(data[[time_name]]==(1:nrow(data)))){
    data[[time_name]]=1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(is.null(offset_name))){
    if (!(offset_name %in% colnames(data))){
      stop("Please make sure that offset_name belongs to colnames(data)")
    }
    if (offset_name %in% covariates | paste0("offset(log(",offset_name,"))") %in% covariates | paste0("offset(",offset_name,")") %in% covariates | paste0("log(",offset_name,")") %in% covariates){
      stop("The offset term should not be included in the formula. Please supply a formula without the offset, and add the name of the offset column using the argument offset_name.")
    }
  }
  if(!is.logical(over_dispersion)){
    stop("over_dispersion must be either TRUE or FALSE")
  }
  if (!(impact_model %in% c("full","level","slope"))){
    stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  }

  # if (offset_name %in% covariates | paste0("offset(log(",offset_name,"))") %in% covariates | paste0("offset(",offset_name,")") %in% covariates | paste0("log(",offset_name,")") %in% covariates){
  #   stop("The offset term should not be included in the formula. Please supply a formula without the offset, and add the name of the offset column using the argument offset_name.")
  # }
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  }
  if (!all(covariates %in% colnames(data))){
    stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  data$indicator <- 0
  data$indicator[intervention_start_ind:n] <- 1
  data$shifted_time <- data[[time_name]]-intervention_start_ind
  # update formula
  if (impact_model=="full"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator + indicator:shifted_time")))
    }
  } else if (impact_model=="level"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator")))
    }
  } else if (impact_model=="slope"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name,"+ indicator:shifted_time")))
    }
  }
  if (!is.null(offset_name)){
    offset_name <- sym(offset_name)
    form_update <-update(form_update, as.formula(paste0("~ . + offset(log(",offset_name,"))")))
  }
  if(!is.logical(keep_significant_fourier)){
    stop("keep_significant_fourier must be either TRUE or FALSE")
  }
  if (!(freq>1 & freq<=n)){
    stop("Please make sure that the freq is a value greater than 1 and less or equal to nrow(data)")
  }
  t_star <- which(data$indicator==1, arr.ind = TRUE)[1]
  #estimate Fourier terms based only on pre-intervention data
  ts_wo_intervention <- ts(data[[response]][1:(t_star-1)],frequency = freq)
  fourier_mat <- fourier(ts_wo_intervention,K=freq/2)
  fourier_mat <- rbind(fourier_mat,fourier_mat[(t_star%%freq):(n-t_star+t_star%%freq),])
  colnames(fourier_mat) <- gsub("-",".", colnames(fourier_mat))
  fourier_df <- data.frame(fourier_mat)
  colnames(fourier_df) <- colnames(fourier_mat)
  data <- cbind(data,fourier_df)
  form_update_full <-update(form_update, as.formula(paste0("~ . +", paste(colnames(fourier_mat), collapse= "+"))))
  # fit model
  if(!isTRUE(over_dispersion)){
    set.seed(1)
    model <- glm(form_update_full,data,family=poisson)
  } else {
    set.seed(1)
    model <- glm(form_update_full,data,family=quasipoisson)
  }
  if(isTRUE(keep_significant_fourier)){
    all_significant_terms <- summary(model)$coeff[,4] < 0.05
    all_fourier_terms <- all_significant_terms[!names(all_significant_terms) %in% c("(Intercept)","time","indicator","indicator:shifted_time")]
    selected_fourier <- names(all_fourier_terms)[all_fourier_terms == TRUE]
    form_update_sig <-update(form_update, as.formula(paste0("~ . +", paste(selected_fourier, collapse= "+"))))
    if(!isTRUE(over_dispersion)){
      set.seed(1)
      model <- glm(form_update_sig,data,family=poisson)
    } else {
      set.seed(1)
      model <- glm(form_update_sig,data,family=quasipoisson)
    }
  }
  b <- (n-t_star)/2
  beta_intervention <- model$coefficients["indicator"]
  beta_interaction <- model$coefficients["indicator:shifted_time"]
  RR <- exp(beta_intervention+beta_interaction*b)
  cov_mat <- vcov(model)
  if (impact_model=="full"){
    beta_intervention <- model$coefficients["indicator"]
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    RR <- exp(beta_intervention+beta_interaction*b)
    var_lin <- cov_mat["indicator","indicator"]+b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]+b*2*cov_mat["indicator","indicator:shifted_time"]
    lin_est <- beta_intervention+b*beta_interaction
  } else if (impact_model=="level"){
    beta_intervention <- model$coefficients["indicator"]
    RR <- exp(beta_intervention)
    var_lin <- cov_mat["indicator","indicator"]
    lin_est <- beta_intervention
  }
  else if (impact_model=="slope"){
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    RR <- exp(beta_interaction*b)
    var_lin <- b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]#+b*2*cov_mat["indicator","indicator:shifted_time"]
    lin_est <- b*beta_interaction
  }
  CI_lin <- c(lin_est-1.96*sqrt(var_lin),lin_est+1.96*sqrt(var_lin))
  CI_RR <- exp(CI_lin)
  TS <- lin_est/sqrt(var_lin)
  p_value <- round(2*pnorm(abs(TS),lower.tail = FALSE),2)
  ret_vec <-c(RR,CI_RR,p_value)
  names(ret_vec) <- c("RR", "2.5% CI", "97.5% CI","P-value")
  s <- summary(model)
  print(s)
  print(ret_vec)
  s[["RR"]] <- ret_vec #or model[["RR"]] <- ret_vec and then return(model)

  #add predictions
  predA <- predict(model,type="response")
  data$pred <- predA# *100/data[[offset_name]]
  if (isTRUE(counterfactual)){
    # Compute counterfactual values
    new_data <- data
    new_data$indicator <- 0
    pred_counter <- predict(model,type="response",newdata = new_data)
    data$predC <- pred_counter#*100/data[[offset_name]]
    data$predC[data$indicator!=1] <- NA
  }
  data <- as_tibble(data)

  # #draw
  # Group <- levels( factor(c("Counterfactual","Fitted values")))
  # response <- sym(response)
  # p2 <- ggplot(data = data , aes(y = !!response*100/data[[offset_name]], x = dt)) +
  #   geom_point()  + ylab(y_lab) + xlab("Date") +
  #   geom_line(aes(y = predC, x = dt, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth)+
  #   geom_line(aes(y = pred, x = dt, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
  #   scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
  #   scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
  #   theme(legend.key=element_blank())
  # p2
  return(list(model=model, model_summary=s, data=data))
}

#' ITS analysis for continuous outcomes with no seasonality
#' \code{its_lm_wo_seas} fits a linear regression model to an ITS, and returns the model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_lm_wo_seas(data=data,form=form,time_name = "time",intervention_start_ind=intervention_start_ind, impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_lm_wo_seas <- function(data,form,time_name, intervention_start_ind,impact_model="full", counterfactual="FALSE"){
  pred <- predC <- NULL
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if(!(is_formula(form))){
    stop("Please make sure that form is a formula object")
  }
  if (!(time_name %in% colnames(data))){
    stop("Please make sure that time_name belongs to colnames(data)")
  }
  if(!all(data[[time_name]]==(1:nrow(data)))){
    data[[time_name]]=1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))){
    stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  }
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  }
  if (!all(covariates %in% colnames(data))){
    stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  data$indicator <- 0
  data$indicator[intervention_start_ind:n] <- 1
  data$shifted_time <- data[[time_name]]-intervention_start_ind
  # update formula
  if (impact_model=="full"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator + indicator:shifted_time")))
    }
  } else if (impact_model=="level"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator")))
    }
  } else if (impact_model=="slope"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name,"+ indicator:shifted_time")))
    }
  }

  # fit model
  model <- lm(form_update,data)
  cov_mat <- vcov(model)

  #Mean difference
  t_star <- which(data$indicator==1, arr.ind = TRUE)[1]
  b <- (n-t_star)/2
  if (impact_model=="full"){
    beta_intervention <- model$coefficients["indicator"]
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    mean_diff <- beta_intervention+beta_interaction*b
    var_lin <- cov_mat["indicator","indicator"]+b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]+b*2*cov_mat["indicator","indicator:shifted_time"]
  } else if (impact_model=="level"){
    beta_intervention <- model$coefficients["indicator"]
    mean_diff <- beta_intervention
    var_lin <- cov_mat["indicator","indicator"]
  } else if (impact_model=="slope"){
    beta_interaction <- model$coefficients["indicator:shifted_time"]
    mean_diff <- beta_interaction*b
    var_lin <- b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]#+b*2*cov_mat["indicator","indicator:shifted_time"]
  }

  CI <- c(mean_diff-1.96*sqrt(var_lin),mean_diff+1.96*sqrt(var_lin))
  TS <- mean_diff/sqrt(var_lin)
  p_value <- round(2*pnorm(abs(TS),lower.tail = FALSE),2)
  ret_vec <-c(mean_diff, CI ,p_value)
  names(ret_vec) <- c("Mean difference", "2.5% CI", "97.5% CI","P-value")
  summar <- summary(model)
  # print(summar)
  # print(ret_vec)
  summar[["Mean difference"]] <- ret_vec #or model[["RR"]] <- ret_vec and then return(model)

  #add predictions and compute effect size
  predA <- predict(model,type="response")
  data$pred <- predA
  intervention_predictions <- data$pred[data$indicator==1]
  s1<- var(intervention_predictions)
  # Compute counterfactual values
  new_data <- data
  new_data$indicator <- 0
  pred_counter <- predict(model,type="response",newdata = new_data)
  data$predC <- pred_counter
  data$predC[data$indicator!=1] <- NA
  data <- as_tibble(data)
  counterfacual_predictions <- na.omit(data$predC)
  s2 <- var(counterfacual_predictions)
  s <- sqrt((s1+s2)/2) #pooled
  # s <- sqrt(var(data[[response]][data$indicator==1])) #SD of observed values for the intervention period
  # s <- sqrt(var(intervention_predictions-counterfacual_predictions)) #SD of difference
  # s <- sqrt(s2) #SD of control
  cohen_d <- mean_diff/s

  #bootstrap for CI and p-value
  #CI
  betas <- model$coefficients
  n_boot_samples <- 2000
  cohen_d_mat <- matrix(NA,n_boot_samples,2)
  model_copy <- model
  model_copy_null <- model
  for (i in 1:n_boot_samples){
    sampled_beta <- mvrnorm(n=1, mu=betas, Sigma=cov_mat )
    sampled_beta_null <- sampled_beta
    if (impact_model=="full"){
      sampled_beta_null[c("indicator","indicator:shifted_time")] <- sampled_beta[c("indicator","indicator:shifted_time")]-betas[c("indicator","indicator:shifted_time")]
    } else if (impact_model=="level"){
      sampled_beta_null["indicator"] <- sampled_beta["indicator"]-betas["indicator"]
    } else if (impact_model=="slope"){
      sampled_beta_null["indicator:shifted_time"] <- sampled_beta["indicator:shifted_time"]-betas["indicator:shifted_time"]
    }
    model_copy$coefficients <- sampled_beta
    fitted_values <- predict(model_copy,type="response")[intervention_start_ind:n]
    counterfactuals <-  predict(model_copy,type="response",newdata = new_data)[intervention_start_ind:n]
    var_fitted <- var(fitted_values)
    var_counter <- var(counterfactuals)
    Sp <- sqrt((var_fitted+var_counter)/2)
    if (impact_model=="full"){
      MD <- sampled_beta["indicator"]+b*sampled_beta["indicator:shifted_time"]
    } else if (impact_model=="level"){
      MD <- sampled_beta["indicator"]
    } else if (impact_model=="slope"){
      MD <- b*sampled_beta["indicator:shifted_time"]
    }
    d <- MD/Sp
    cohen_d_mat[i,1] <- d
    model_copy$coefficients <- sampled_beta_null
    fitted_values <- predict(model_copy,type="response")[intervention_start_ind:n]
    counterfactuals <-  predict(model_copy,type="response",newdata = new_data)[intervention_start_ind:n]
    var_fitted <- var(fitted_values)
    var_counter <- var(counterfactuals)
    Sp <- sqrt((var_fitted+var_counter)/2)
    if (impact_model=="full"){
      MD <- sampled_beta_null["indicator"]+b*sampled_beta_null["indicator:shifted_time"]
    } else if (impact_model=="level"){
      MD <- sampled_beta_null["indicator"]
    } else if (impact_model=="slope"){
      MD <- b*sampled_beta_null["indicator:shifted_time"]
    }
    d <- MD/Sp
    cohen_d_mat[i,2] <- d
  }
  CI_low <- quantile(cohen_d_mat[,1],0.025)
  CI_up <- quantile(cohen_d_mat[,1],0.975)
  p_value_d <- mean(abs(cohen_d_mat[,2])>abs(cohen_d))

  # #P-value
  # null_betas <- betas
  # if (impact_model=="full"){
  #   null_betas[c("indicator","indicator:shifted_time")] <- c(0,0)
  # } else if (impact_model=="level"){
  #   null_betas["indicator"] <- 0
  # } else if (impact_model=="slope"){
  #   null_betas["indicator:shifted_time"] <- 0
  # }
  # cohen_d_vector <- matrix(NA,n_boot_samples,1)
  # for (i in 1:n_boot_samples){
  #   sampled_null_beta <- mvrnorm(n=1, mu=null_betas, Sigma=cov_mat )
  #   model_copy$coefficients <- sampled_null_beta
  #   fitted_values <- predict(model_copy,type="response")[intervention_start_ind:n]
  #   counterfactuals <-  predict(model_copy,type="response",newdata = new_data)[intervention_start_ind:n]
  #   var_fitted <- var(fitted_values)
  #   var_counter <- var(counterfactuals)
  #   Sp <- sqrt((var_fitted+var_counter)/2)
  #   if (impact_model=="full"){
  #     MD <- sampled_null_beta["indicator"]+b*sampled_null_beta["indicator:shifted_time"]
  #   } else if (impact_model=="level"){
  #     MD <- sampled_null_beta["indicator"]
  #   }
  #   else if (impact_model=="slope"){
  #     MD <- b*sampled_null_beta["indicator:shifted_time"]
  #   }
  #   d <- MD/Sp
  #   cohen_d_vector[i] <- d
  # }
  # p_value_d <- mean(abs(cohen_d_vector)>abs(cohen_d))
  ret_vec2 <-c(cohen_d, CI_low, CI_up ,p_value_d)
  names(ret_vec2) <- c("Cohen's d", "2.5% CI", "97.5% CI","P-value")
  print(summar)
  print(ret_vec)
  print(ret_vec2)
  summar[["Cohen's d"]] <- ret_vec2

  if (!(isTRUE(counterfactual))){
    data$predC <- NULL
  }
  # N <- 2*length(intervention_predictions)-2
  # J <- gamma(N/2) / (sqrt(N/2) * gamma((N-1)/2))
  # bias_corrected_d <- J*cohen_d
  # sig_d <- sqrt((8+bias_corrected_d^2)/(4*length(intervention_predictions)))
  # lower <- bias_corrected_d-1.96*sig_d
  # upper <- bias_corrected_d+1.96*sig_d
  # TS_d <- bias_corrected_d/sig_d
  # p_value_d <- round(2*pnorm(abs(TS_d),lower.tail = FALSE),2)



  # #draw
  # Group <- levels( factor(c("Counterfactual","Fitted values")))
  # response <- sym(response)
  # p2 <- ggplot(data = data , aes(y = !!response, x = dt)) +
  #   geom_point()  + ylab(y_lab) + xlab("Date") +
  #   geom_line(aes(y = predC, x = dt, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth)+
  #   geom_line(aes(y = pred, x = dt, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
  #   scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
  #   scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
  #   theme(legend.key=element_blank())
  # p2
  return(list(model=model, model_summary=summar, data=data))
}

#' ITS analysis for continuous outcomes with seasonal adjustments via Fourier terms.
#'
#' \code{its_lm_fourier} fits a linear regression model asjusted to seasonality, and returns the model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param freq A positive integer describing the frequency of the time series.
#' @param keep_significant_fourier Logical - indicating whether only the significant Fourier terms should be considered. Default is TRUE and then the model is fitted twice; once to obtain the significant Fourier terms, and second time keeping only the significant Fourier terms. If FALSE, then all the Fourier terms are used.
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_lm_fourier(data=data,form=form,time_name = "time",intervention_start_ind=intervention_start_ind,freq=12, keep_significant_fourier=TRUE, impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_lm_fourier <- function(data, form, time_name, intervention_start_ind,freq, keep_significant_fourier=TRUE,impact_model="full", counterfactual="FALSE"){
  pred <- predC <- NULL
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if(!(is_formula(form))){
    stop("Please make sure that form is a formula object")
  }
  if (!(time_name %in% colnames(data))){
    stop("Please make sure that time_name belongs to colnames(data)")
  }
  if(!all(data[[time_name]]==(1:nrow(data)))){
    data[[time_name]]=1:nrow(data)
    print("time_name column was overwritten with the values 1:nrow(data)")
  }
  if (!(impact_model %in% c("full","level","slope"))){
    stop("Please enter a valid impact_model name. Possible impact models are \"full\", \"level\", and \"slope\".")
  }
  vars <- all.vars(form)
  response <- vars[1]
  covariates <- vars[-1]
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable on the left hand side of the formula object belongs to colnames(data)")
  }
  if (!all(covariates %in% colnames(data))){
    stop("Please make sure that the covariates on the right hand side of the formula object belong to colnames(data)")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  data$indicator <- 0
  data$indicator[intervention_start_ind:n] <- 1
  data$shifted_time <- data[[time_name]]-intervention_start_ind
  # update formula
  if (impact_model=="full"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator + indicator:shifted_time")))
    }
  } else if (impact_model=="level"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name," + indicator")))
    }
  } else if (impact_model=="slope"){
    if (time_name %in% covariates){
      form_update <-update(form,    as.formula(paste0("~ . + indicator:shifted_time")))
    } else {
      form_update <-update(form,    as.formula(paste0("~ . +", time_name,"+ indicator:shifted_time")))
    }
  }
  if(!is.logical(keep_significant_fourier)){
    stop("keep_significant_fourier must be either TRUE or FALSE")
  }
  if (!(freq>1 & freq<=n)){
    stop("Please make sure that the freq is a value greater than 1 and less or equal to nrow(data)")
  }
  t_star <- which(data$indicator==1, arr.ind = TRUE)[1]
  #estimate Fourier terms based only on pre-intervention data
  ts_wo_intervention <- ts(data[[response]][1:(t_star-1)],frequency = freq)
  fourier_mat <- fourier(ts_wo_intervention,K=freq/2)
  fourier_mat <- rbind(fourier_mat,fourier_mat[(t_star%%freq):(n-t_star+t_star%%freq),])
  colnames(fourier_mat) <- gsub("-",".", colnames(fourier_mat))
  fourier_df <- data.frame(fourier_mat)
  colnames(fourier_df) <- colnames(fourier_mat)
  data <- cbind(data,fourier_df)
  form_update_full <-update(form_update, as.formula(paste0("~ . +", paste(colnames(fourier_mat), collapse= "+"))))
  # fit model
  model <- lm(form_update_full,data)
  if(isTRUE(keep_significant_fourier)){
    all_significant_terms <- summary(model)$coeff[,4] < 0.05
    all_fourier_terms <- all_significant_terms[!names(all_significant_terms) %in% c("(Intercept)","time","indicator","indicator:shifted_time")]
    selected_fourier <- names(all_fourier_terms)[all_fourier_terms == TRUE]
    form_update_sig <-update(form_update, as.formula(paste0("~ . +", paste(selected_fourier, collapse= "+"))))
    model <- lm(form_update_sig,data)
  }
  b <- (n-t_star)/2
  beta_intervention <- model$coefficients["indicator"]
  beta_interaction <- model$coefficients["indicator:shifted_time"]
  mean_diff <- beta_intervention+beta_interaction*b
  cov_mat <- vcov(model)
  var_lin <- cov_mat["indicator","indicator"]+b^2*cov_mat["indicator:shifted_time","indicator:shifted_time"]+b*2*cov_mat["indicator","indicator:shifted_time"]
  #lin_est <- beta_intervention+b*beta_interaction
  CI <- c(mean_diff-1.96*sqrt(var_lin),mean_diff+1.96*sqrt(var_lin))
  TS <- mean_diff/sqrt(var_lin)
  p_value <- round(2*pnorm(abs(TS),lower.tail = FALSE),2)
  ret_vec <-c(mean_diff, CI ,p_value)
  names(ret_vec) <- c("Mean difference", "2.5% CI", "97.5% CI","P-value")
  summar <- summary(model)
  # print(summar)
  # print(ret_vec)
  summar[["Mean difference"]] <- ret_vec #or model[["RR"]] <- ret_vec and then return(model)

  #add predictions and compute effect size
  predA <- predict(model,type="response")
  data$pred <- predA
  intervention_predictions <- data$pred[data$indicator==1]
  s1<- var(intervention_predictions)
  # Compute counterfactual values
  new_data <- data
  new_data$indicator <- 0
  pred_counter <- predict(model,type="response",newdata = new_data)
  data$predC <- pred_counter
  data$predC[data$indicator!=1] <- NA
  data <- as_tibble(data)
  counterfacual_predictions <- na.omit(data$predC)
  s2 <- var(counterfacual_predictions)
  s <- sqrt((s1+s2)/2)
  # s <- sqrt(var(data[[response]][data$indicator==1]))
  # s <- sqrt(var(intervention_predictions-counterfacual_predictions))
  # s <- sqrt(s2)
  cohen_d <- mean_diff/s

  #bootstrap for CI and p-value
  #CI
  betas <- model$coefficients
  n_boot_samples <- 2000
  cohen_d_mat <- matrix(NA,n_boot_samples,2)
  model_copy <- model
  model_copy_null <- model
  for (i in 1:n_boot_samples){
    sampled_beta <- mvrnorm(n=1, mu=betas, Sigma=cov_mat )
    sampled_beta_null <- sampled_beta
    if (impact_model=="full"){
      sampled_beta_null[c("indicator","indicator:shifted_time")] <- sampled_beta[c("indicator","indicator:shifted_time")]-betas[c("indicator","indicator:shifted_time")]
    } else if (impact_model=="level"){
      sampled_beta_null["indicator"] <- sampled_beta["indicator"]-betas["indicator"]
    } else if (impact_model=="slope"){
      sampled_beta_null["indicator:shifted_time"] <- sampled_beta["indicator:shifted_time"]-betas["indicator:shifted_time"]
    }
    model_copy$coefficients <- sampled_beta
    fitted_values <- predict(model_copy,type="response")[intervention_start_ind:n]
    counterfactuals <-  predict(model_copy,type="response",newdata = new_data)[intervention_start_ind:n]
    var_fitted <- var(fitted_values)
    var_counter <- var(counterfactuals)
    Sp <- sqrt((var_fitted+var_counter)/2)
    if (impact_model=="full"){
      MD <- sampled_beta["indicator"]+b*sampled_beta["indicator:shifted_time"]
    } else if (impact_model=="level"){
      MD <- sampled_beta["indicator"]
    } else if (impact_model=="slope"){
      MD <- b*sampled_beta["indicator:shifted_time"]
    }
    d <- MD/Sp
    cohen_d_mat[i,1] <- d
    model_copy$coefficients <- sampled_beta_null
    fitted_values <- predict(model_copy,type="response")[intervention_start_ind:n]
    counterfactuals <-  predict(model_copy,type="response",newdata = new_data)[intervention_start_ind:n]
    var_fitted <- var(fitted_values)
    var_counter <- var(counterfactuals)
    Sp <- sqrt((var_fitted+var_counter)/2)
    if (impact_model=="full"){
      MD <- sampled_beta_null["indicator"]+b*sampled_beta_null["indicator:shifted_time"]
    } else if (impact_model=="level"){
      MD <- sampled_beta_null["indicator"]
    } else if (impact_model=="slope"){
      MD <- b*sampled_beta_null["indicator:shifted_time"]
    }
    d <- MD/Sp
    cohen_d_mat[i,2] <- d
  }
  CI_low <- quantile(cohen_d_mat[,1],0.025)
  CI_up <- quantile(cohen_d_mat[,1],0.975)
  p_value_d <- mean(abs(cohen_d_mat[,2])>abs(cohen_d))

  ret_vec2 <-c(cohen_d, CI_low, CI_up ,p_value_d)
  names(ret_vec2) <- c("Cohen's d", "2.5% CI", "97.5% CI","P-value")
  print(summar)
  print(ret_vec)
  print(ret_vec2)
  summar[["Cohen's d"]] <- ret_vec2

  # N <- 2*length(intervention_predictions)-2
  # J <- gamma(N/2) / (sqrt(N/2) * gamma((N-1)/2))
  # bias_corrected_d <- J*cohen_d
  # sig_d <- sqrt((8+bias_corrected_d^2)/(4*length(intervention_predictions)))
  # lower <- bias_corrected_d-1.96*sig_d
  # upper <- bias_corrected_d+1.96*sig_d
  # TS_d <- bias_corrected_d/sig_d
  # p_value_d <- round(2*pnorm(abs(TS_d),lower.tail = FALSE),2)
  #
  # ret_vec2 <-c(bias_corrected_d, lower, upper ,p_value_d)
  # names(ret_vec2) <- c("Adjusted Cohen's d", "2.5% CI", "97.5% CI","P-value")
  # print(summar)
  # print(ret_vec)
  # print(ret_vec2)
  # summar[["Cohen's d"]] <- ret_vec2

  if (!(isTRUE(counterfactual))){
    data$predC <- NULL
  }
  # #draw
  # Group <- levels( factor(c("Counterfactual","Fitted values")))
  # response <- sym(response)
  # p2 <- ggplot(data = data , aes(y = !!response, x = dt)) +
  #   geom_point()  + ylab(y_lab) + xlab("Date") +
  #   geom_line(aes(y = predC, x = dt, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth)+
  #   geom_line(aes(y = pred, x = dt, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
  #   scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
  #   scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
  #   theme(legend.key=element_blank())
  # p2
  return(list(model=model, model_summary=summar, data=data))
}


#' ITS analysis for continuous outcomes
#'
#' \code{its_lm} fits a linear regression model to an ITS, and returns the model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param freq A positive integer describing the frequency of the time series.
#' @param seasonality A string specifying whether seasonality should be considered. Possible options include "none" corresponding to no seasonal adjustment, "full" corresponding to using freq-1 Fourier terms to model the seasonal component, and "significant" indicating whether only the significant Fourier terms should be considered in the seasonal adjustment. Default value is "none".
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted linear regression model, the summary of the model (including the mean difference and Cohen's d), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_lm(data=data,form=form,time_name = "time",intervention_start_ind=intervention_start_ind,freq=12, seasonality= "none", impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_lm <- function(data, form, time_name, intervention_start_ind,freq, seasonality= "none",impact_model="full", counterfactual=FALSE){
  pred <- predC <- NULL
  if (seasonality=="none"){
    out <- its_lm_wo_seas(data, form, time_name, intervention_start_ind,impact_model, counterfactual)
  } else if (seasonality=="full"){
    keep_significant_fourier=FALSE
    out <- its_lm_fourier(data, form, time_name, intervention_start_ind,freq,keep_significant_fourier, impact_model, counterfactual)
  } else if (seasonality=="significant"){
    keep_significant_fourier=TRUE
    out <- its_lm_fourier(data, form, time_name, intervention_start_ind,freq,keep_significant_fourier, impact_model, counterfactual)
  } else {
    stop("Please enter a valid seasonality argument. Possible seasonality arguments are \"none\", \"full\", and \"significant\".")
  }
  return(out)
}

#' ITS analysis for count outcomes
#'
#' \code{its_poisson} fits a Poisson regression model to an ITS, and returns the model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param form A formula with the response on the left, followed by the ~ operator, and the covariates on the right, separated by + operators. The formula should not contain an offest term.
#' @param offset_name either a string indicating the name of the offset column in the data, or NULL. Default value is NULL
#' @param time_name A string giving the name of the time variable. The time variable may or may not be supplied as a covariate in the formula
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param over_dispersion Logical - indicating whether a quasi-Poisson model should be used to account for over dispersion (TRUE), or not (FLASE). Default value is FALSE.
#' @param freq A positive integer describing the frequency of the time series.
#' @param seasonality A string specifying whether seasonality should be considered. Possible options include "none" corresponding to no seasonal adjustment, "full" corresponding to using freq-1 Fourier terms to model the seasonal component, and "significant" indicating whether only the significant Fourier terms should be considered in the seasonal adjustment. Default value is "none".
#' @param impact_model A string specifying the assumed impact model. Possible options include "full" corresponding to a model including both a level change and a slope change, "level" corresponding to a model including just a level change, and "slope" corresponding to a model including just a slope change. Default value is "full".
#' @param counterfactual Logical - indicating whether the model-based counterfactual values should also be returned as an additional column in the data. Default value is FALSE, in which case the counterfactual values are not returned.
#' @return The function returns a list with three elements: the fitted Poisson regression model, the summary of the model (including the relative risk), and the original data together with the model predictions.
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_poisson(data=data,form=form,offset_name="labour", time_name = "time",intervention_start_ind=intervention_start_ind,over_dispersion=TRUE, freq=12, seasonality= "none", impact_model = "full",counterfactual = TRUE)
#' @importFrom tibble is_tibble as_tibble
#' @importFrom rlang is_formula sym
#' @importFrom stats as.formula glm lm pnorm poisson quasipoisson ts vcov update predict var na.omit quantile
#' @importFrom forecast fourier
#' @importFrom MASS mvrnorm
#' @export
its_poisson <- function(data, form, offset_name=NULL,time_name, intervention_start_ind,over_dispersion=FALSE, freq, seasonality= "none",impact_model="full", counterfactual=FALSE){
  pred <- predC <- NULL
  if (seasonality=="none"){
    out <- its_poisson_wo_seas(data, form, offset_name,time_name, intervention_start_ind,over_dispersion, impact_model, counterfactual)
  } else if (seasonality=="full"){
    keep_significant_fourier=FALSE
    out <- its_poisson_fourier(data, form, offset_name,time_name, intervention_start_ind,over_dispersion, freq, keep_significant_fourier,impact_model, counterfactual)
  } else if (seasonality=="significant"){
    keep_significant_fourier=TRUE
    out <- its_poisson_fourier(data, form, offset_name,time_name, intervention_start_ind,over_dispersion, freq, keep_significant_fourier,impact_model, counterfactual)
  } else {
    stop("Please enter a valid seasonality argument. Possible seasonality arguments are \"none\", \"full\", and \"significant\".")
  }
  return(out)
}

#' Plot ITS fitted values
#'
#' \code{plot_its_lm} uses ggplot2 to plot the model-based fitted values, together with a scatterplot of the observed time series
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param y_lab a string with the y axis label for the predictions
#' @param response The name of the response variable to be plotted in the scatterplot
#' @param date_name A string giving the name of the date column. The date column must be a Date object.
#' @return a ggplot object including a scatterplot of the time series, the predictions line, and the counterfactual predictions in red (if available.
#' @examples
#' data <- unemployed
#' form <- as.formula("percent ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_lm(data=data,form=form, time_name = "time",intervention_start_ind=intervention_start_ind, freq=12, seasonality= "none", impact_model = "full",counterfactual = TRUE)
#' new_data <- fit$data
#' y_lab <- "Unemployment percent"
#' response <- "percent"
#' p <- plot_its_lm(new_data,intervention_start_ind,y_lab,response, "dt")
#' @importFrom lubridate is.Date
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab geom_line draw_key_smooth geom_segment scale_color_manual scale_linetype_manual theme element_blank theme_bw
#' @export
plot_its_lm <- function(data,intervention_start_ind,y_lab,response,date_name){
  pred <- predC <- NULL
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable belongs to colnames(data)")
  }
  if (!(date_name %in% colnames(data))){
    stop("Please make sure that the date_name column belongs to colnames(data)")
  } else if (!(is.Date(data[[date_name]]))){
    stop("Please make sure that the date_name column is a Date object")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  if (!("pred" %in% colnames(data))){
    stop("The data needs to include the column of predictions. Please use as input the data output of the function its_lm()")
  }
  if (!(is.character(y_lab))){
    y_lab <- ""
    print("y_lab is not a string and thus will not by used")
  }
  response <- sym(response)
  xdot <- data[[date_name]][intervention_start_ind]
  date_name <- sym(date_name)
  data_pre <- data[1:(intervention_start_ind-1),]
  data_post <- data[intervention_start_ind:nrow(data),]
  if ("predC" %in% colnames(data)){
    Group <- levels( factor(c("Counterfactual","Fitted values")))
    p <- ggplot(data = data , aes(y = !!response, x = !!date_name)) +
      geom_point()  + ylab(y_lab) + xlab("Date") +
      geom_line(aes(y = predC, x = !!date_name, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth,na.rm=TRUE)+
      geom_line(data = data_pre, aes(y = pred, x = !!date_name, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
      geom_line(data = data_post, aes(y = pred, x = !!date_name, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
      geom_segment(data = data, aes(x=xdot,xend=xdot,y=pred[intervention_start_ind-1],yend=pred[intervention_start_ind]), linetype=3,size=1) +
      scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
      scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
      theme(legend.key=element_blank()) +theme_bw()
  } else {
    p <- ggplot(data = data , aes(y = !!response, x = !!date_name)) +
      geom_point()  + ylab(y_lab) + xlab("Date") +
      #geom_line(aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth,na.rm=TRUE)+
      geom_line(data = data_pre, aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth)+
      geom_line(data = data_post, aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth)+
      geom_segment(data = data, aes(x=xdot,xend=xdot,y=pred[intervention_start_ind-1],yend=pred[intervention_start_ind]), linetype=3,size=1) +
      theme(legend.key=element_blank()) +theme_bw()
  }
  return(p)
}

#' Plot ITS fitted values
#'
#' \code{plot_its_poisson} uses ggplot2 to plot the model-based fitted values, together with a scatterplot of the observed time series
#'
#' @param data The data frame corresponding to the supplied formula, existing of at least 2 variables: (1) the count outcome, and (2) a vector of time points
#' @param intervention_start_ind Numeric - a number between 1 and nrow(data)-1 stating the time point of the start of the intervention
#' @param y_lab a string with the y axis label for the predictions
#' @param response The name of the response variable to be plotted in the scatterplot
#' @param offset_name either a string indicating the name of the offset column used in the ITS Poisson model, or NULL. Default value is NULL. When offfset_name exists, the function plots the response rate, instead of the count itself.
#' @param date_name A string giving the name of the date column. The date column must be a Date object.
#' @return a ggplot object including a scatterplot of the time series, the predictions line, and the counterfactual predictions in red (if available).
#' @examples
#' data <- unemployed
#' form <- as.formula("unemployed ~ time")
#' intervention_start_ind <- which(data$year==2020 & data$month>2| data$year==2021)[1]
#' fit <- its_poisson(data=data,form=form,offset_name="labour", time_name = "time",intervention_start_ind=intervention_start_ind,over_dispersion=TRUE, freq=12, seasonality= "none", impact_model = "full",counterfactual = TRUE)
#' new_data <- fit$data
#' y_lab <- "Unemployment percent"
#' response <- "percent"
#' p <- plot_its_poisson(new_data,intervention_start_ind,y_lab,response, offset_name = "labour", "dt")
#' @importFrom lubridate is.Date
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab geom_line draw_key_smooth geom_segment scale_color_manual scale_linetype_manual theme element_blank theme_bw
#' @export
plot_its_poisson <- function(data,intervention_start_ind,y_lab,response, offset_name = NULL, date_name){
  pred <- predC <- NULL
  if (!(is.data.frame(data) | is_tibble(data))){
    stop("Please make sure that data is either a data frame or a tibble")
  }
  if (!(response %in% colnames(data))){
    stop("Please make sure that the response variable belongs to colnames(data)")
  }
  if (!(date_name %in% colnames(data))){
    stop("Please make sure that the date_name column belongs to colnames(data)")
  } else if (!(is.Date(data[[date_name]]))){
    stop("Please make sure that the date_name column is a Date object")
  }
  n <- nrow(data)
  if (!(intervention_start_ind>1 & intervention_start_ind<=n)){
    stop("Please make sure that intervention_start_ind is a value greater than 1 and less or equal to nrow(data). intervention_start_ind is the index ")
  }
  if (!("pred" %in% colnames(data))){
    stop("The data needs to include the column of predictions. Please use as input the data output of the function its_lm()")
  }
  if (!(is.character(y_lab))){
    y_lab <- ""
    print("y_lab is not a string and thus will not by used")
  }
  if (!(is.null(offset_name))){
    if (!(offset_name %in% colnames(data))){
      stop("Please make sure that offset_name belongs to colnames(data)")
    } else {
      data$pred <- data$pred*100/data[[offset_name]]
      data[[response]] <- data[[response]]*100/data[[offset_name]]
      if ("predC" %in% colnames(data)){
        data$predC <- data$predC*100/data[[offset_name]]
      }
    }
  }
  response <- sym(response)
  xdot <- data[[date_name]][intervention_start_ind]
  date_name <- sym(date_name)
  data_pre <- data[1:(intervention_start_ind-1),]
  data_post <- data[intervention_start_ind:nrow(data),]
  if ("predC" %in% colnames(data)){
    Group <- levels( factor(c("Counterfactual","Fitted values")))
    p <- ggplot(data = data , aes(y = !!response, x = !!date_name)) +
      geom_point()  + ylab(y_lab) + xlab("Date") +
      geom_line(aes(y = predC, x = !!date_name, color="Counterfactual",linetype = "Counterfactual"), size=1,key_glyph = draw_key_smooth,na.rm=TRUE)+
      geom_line(data = data_pre, aes(y = pred, x = !!date_name, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
      geom_line(data = data_post, aes(y = pred, x = !!date_name, color="Fitted values",linetype = "Fitted values"), size=1,key_glyph = draw_key_smooth)+
      geom_segment(data = data, aes(x=xdot,xend=xdot,y=pred[intervention_start_ind-1],yend=pred[intervention_start_ind]), linetype=3,size=1) +
      scale_color_manual("", values = c("Counterfactual"="red","Fitted values"="black"),breaks = Group) +
      scale_linetype_manual("",values = c("Counterfactual"=2,"Fitted values"=1),breaks = Group) +
      theme(legend.key=element_blank()) +theme_bw()
  } else {
    p <- ggplot(data = data , aes(y = !!response, x = !!date_name)) +
      geom_point()  + ylab(y_lab) + xlab("Date") +
      #geom_line(aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth,na.rm=TRUE)+
      geom_line(data = data_pre, aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth)+
      geom_line(data = data_post, aes(y = pred, x = !!date_name), size=1,key_glyph = draw_key_smooth)+
      geom_segment(data = data, aes(x=xdot,xend=xdot,y=pred[intervention_start_ind-1],yend=pred[intervention_start_ind]), linetype=3,size=1) +
      theme(legend.key=element_blank()) +theme_bw()
  }
  return(p)
}
