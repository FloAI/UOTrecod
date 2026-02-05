#' OT_unbalanced_joint
#'
#' @param datab data.frame containing ID, Y, Z and covariates
#' @param index_DB_Y_Z Vector of 3 indexes (ID, Y, Z)
#' @param tau Penalty parameter for unbalanced OT (relaxation of marginals)
#' @param lambda.reg Regularization parameter for R-JOINT logic
#' @param dist.choice "E", "M", "G", or "H"
#' @param which.DB "A", "B", or "BOTH"
OT_unbalanced_joint <- function(datab, index_DB_Y_Z = 1:3,
                                nominal = NULL, ordinal = NULL, logic = NULL,
                                convert.num = NULL, convert.class = NULL,
                                dist.choice = "E", percent.knn = 1,
                                tau = 1.0, lambda.reg = 0.0, prox.X = 0.30,
                                solvR = "glpk", which.DB = "BOTH") {

  tstart <- Sys.time()

  # 1. Validation and Pre-processing
  if (!(dist.choice %in% c("M", "E", "G", "H"))) {
    dist.choice <- substr(dist.choice, 1, 1)
  }

  message("---------------------------------------")
  message("OT UNBALANCED JOINT PROCEDURE")
  message("Distance: ", dist.choice, " | Tau: ", tau)
  message("---------------------------------------")

  dataB <- transfo_dist(datab, index_DB_Y_Z = index_DB_Y_Z, quanti = convert.num,
                        nominal = nominal, ordinal = ordinal, logic = logic,
                        convert_num = convert.num, convert_class = convert.class,
                        prep_choice = dist.choice)

  inst <- proxim_dist(dataB, norm = dist.choice, prox = 0)

  # 2. Extract unique covariate profiles and handle potential empty results
  Xvalues <- unique(inst$Xobserv)
  if (is.null(Xvalues) || nrow(Xvalues) == 0) {
    stop("Error: No covariate profiles (X) were found. Check your index_DB_Y_Z and covariate column indexes.")
  }

  nA <- inst$nA; nB <- inst$nB; Y <- inst$Y; Z <- inst$Z
  indXA <- inst$indXA; indXB <- inst$indXB; nbX <- length(indXA)
  Xobserv <- inst$Xobserv; Zobserv <- inst$Zobserv; Yobserv <- inst$Yobserv
  ID_prof <- paste0("P_", 1:nrow(Xvalues))

  # 3. Distance Matrix Calculation for Regularization
  if (dist.choice == "G") {
    dist_X <- StatMatch::gower.dist(Xvalues, Xvalues)
  } else {
    # Coerce to matrix to ensure proxy::dist doesn't return an empty object
    dist_X <- as.matrix(proxy::dist(as.matrix(Xvalues),
                                    method = switch(dist.choice, M="manhattan", E="euclidean", H="hamming")))
  }

  tol.X <- prox.X * max(dist_X, na.rm = TRUE)
  voisins_X <- dist_X <= tol.X
  ind_voisins <- lapply(1:nrow(voisins_X), function(i) which(voisins_X[i, ]))

  C <- avg_dist_closest(inst, percent_closest = percent.knn)[[1]]

  # 4. Marginal Estimators
  estim_XA <- sapply(indXA, function(x) length(x)/nA)
  estim_XB <- sapply(indXB, function(x) length(x)/nB)

  estim_XA_YA <- lapply(1:nbX, function(x) sapply(Y, function(y) length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA))
  estim_XB_ZB <- lapply(1:nbX, function(x) sapply(Z, function(z) length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB))

  # --- OPTIMIZATION FOR DATABASE A (Predict Z in A) ---
  gammaA_val <- NULL; estimatorZA <- NULL; DATA1_OT <- NULL
  if (which.DB %in% c("A", "BOTH")) {
    modelA <- MIPModel() %>%
      add_variable(gammaA[x, y, z], x = 1:nbX, y = Y, z = Z, type = "continuous", lb = 0) %>%
      add_variable(slack_YA[x, y], x = 1:nbX, y = Y, lb = 0) %>%
      add_variable(slack_ZA[x, z], x = 1:nbX, z = Z, lb = 0) %>%
      add_variable(reg_absA[x1, x2, y, z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z, lb = 0) %>%
      set_objective(
        sum_expr(C[y, z] * gammaA[x, y, z], x = 1:nbX, y = Y, z = Z) +
        tau * (sum_expr(slack_YA[x, y], x = 1:nbX, y = Y) + sum_expr(slack_ZA[x, z], x = 1:nbX, z = Z)) +
        sum_expr((lambda.reg/length(ind_voisins[[x1]])) * reg_absA[x1, x2, y, z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z), "min"
      ) %>%
      add_constraint(sum_expr(gammaA[x, y, z], z = Z) + slack_YA[x, y] >= estim_XA_YA[[x]][y], x = 1:nbX, y = Y) %>%
      add_constraint(sum_expr(gammaA[x, y, z], y = Y) + slack_ZA[x, z] >= (estim_XB_ZB[[x]][z] * estim_XA[[x]]) / max(1e-6, estim_XB[[x]]), x = 1:nbX, z = Z) %>%
      add_constraint(reg_absA[x1, x2, y, z] + gammaA[x2, y, z]/max(1e-6, estim_XA[x2]) >= gammaA[x1, y, z]/max(1e-6, estim_XA[x1]), x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z) %>%
      solve_model(with_ROI(solver = solvR))

    solA <- get_solution(modelA, gammaA[x, y, z])
    gammaA_val <- array(solA$value, dim = c(nbX, length(Y), length(Z)))
    estimatorZA <- array(1/length(Z), dim = c(nbX, length(Y), length(Z)))
    for(x in 1:nbX) {
      for(y in Y) {
        mass <- sum(gammaA_val[x, y, ])
        if(mass > 1e-6) estimatorZA[x, y, ] <- gammaA_val[x, y, ] / mass
      }
    }
    probaZindivA <- matrix(0, nA, length(Z))
    for(x in 1:nbX) for(i in indXA[[x]]) probaZindivA[i,] <- estimatorZA[x, Yobserv[i], ]
    predZA <- apply(probaZindivA, 1, which.max)
    DATA1_OT <- dataB[dataB[, 1] == unique(dataB[, 1])[1], ]
    DATA1_OT$OTpred <- factor(levels(dataB[, 3])[predZA], levels = levels(dataB[, 3]))
  }

  # --- OPTIMIZATION FOR DATABASE B (Predict Y in B) ---
  gammaB_val <- NULL; estimatorYB <- NULL; DATA2_OT <- NULL
  if (which.DB %in% c("B", "BOTH")) {
    modelB <- MIPModel() %>%
      add_variable(gammaB[x, y, z], x = 1:nbX, y = Y, z = Z, type = "continuous", lb = 0) %>%
      add_variable(slack_YB[x, y], x = 1:nbX, y = Y, lb = 0) %>%
      add_variable(slack_ZB[x, z], x = 1:nbX, z = Z, lb = 0) %>%
      add_variable(reg_absB[x1, x2, y, z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z, lb = 0) %>%
      set_objective(
        sum_expr(C[y, z] * gammaB[x, y, z], x = 1:nbX, y = Y, z = Z) +
        tau * (sum_expr(slack_YB[x, y], x = 1:nbX, y = Y) + sum_expr(slack_ZB[x, z], x = 1:nbX, z = Z)) +
        sum_expr((lambda.reg/length(ind_voisins[[x1]])) * reg_absB[x1, x2, y, z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z), "min"
      ) %>%
      add_constraint(sum_expr(gammaB[x, y, z], y = Y) + slack_ZB[x, z] >= estim_XB_ZB[[x]][z], x = 1:nbX, z = Z) %>%
      add_constraint(sum_expr(gammaB[x, y, z], z = Z) + slack_YB[x, y] >= (estim_XA_YA[[x]][y] * estim_XB[[x]]) / max(1e-6, estim_XA[[x]]), x = 1:nbX, y = Y) %>%
      add_constraint(reg_absB[x1, x2, y, z] + gammaB[x2, y, z]/max(1e-6, estim_XB[x2]) >= gammaB[x1, y, z]/max(1e-6, estim_XB[x1]), x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z) %>%
      solve_model(with_ROI(solver = solvR))

    solB <- get_solution(modelB, gammaB[x, y, z])
    gammaB_val <- array(solB$value, dim = c(nbX, length(Y), length(Z)))
    estimatorYB <- array(1/length(Y), dim = c(nbX, length(Z), length(Y)))
    for(x in 1:nbX) {
      for(z in Z) {
        mass <- sum(gammaB_val[x, , z])
        if(mass > 1e-6) estimatorYB[x, z, ] <- gammaB_val[x, , z] / mass
      }
    }
    probaYindivB <- matrix(0, nB, length(Y))
    for(x in 1:nbX) for(i in indXB[[x]]) probaYindivB[i,] <- estimatorYB[x, Zobserv[i + nA], ]
    predYB <- apply(probaYindivB, 1, which.max)
    DATA2_OT <- dataB[dataB[, 1] == unique(dataB[, 1])[2], ]
    DATA2_OT$OTpred <- factor(levels(dataB[, 2])[predYB], levels = levels(dataB[, 2]))
  }

  tend <- Sys.time()

  # Final Object Creation
  res <- list(
    time_exe = difftime(tend, tstart),
    gamma_A = if(!is.null(gammaA_val)) apply(gammaA_val, c(2,3), sum) else NULL,
    gamma_B = if(!is.null(gammaB_val)) apply(gammaB_val, c(2,3), sum) else NULL,
    profile = data.frame(ID = ID_prof, Xvalues),
    res_prox = inst,
    estimatorZA = estimatorZA,
    estimatorYB = estimatorYB,
    DATA1_OT = if(is.null(DATA1_OT)) dataB[dataB[, 1] == unique(dataB[, 1])[1], ] else DATA1_OT,
    DATA2_OT = if(is.null(DATA2_OT)) dataB[dataB[, 1] == unique(dataB[, 1])[2], ] else DATA2_OT
  )
  class(res) <- "otres"
  return(res)
}