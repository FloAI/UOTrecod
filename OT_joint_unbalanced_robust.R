OT_Metric_Robust_KL_joint <- function(datab, index_DB_Y_Z = 1:3,
                                      nominal = NULL, ordinal = NULL, logic = NULL,
                                      convert.num = NULL, convert.class = NULL,
                                      tau_A = 1.0, tau_B = 1.0, epsilon = 0.05,
                                      k_huber = 1.345, which.DB = "BOTH") {

  tstart <- Sys.time()

  # --- 1. PRE-PROCESSING ---
  # Standard otrecod logic to handle factors and variable types
  dataB <- transfo_dist(datab, index_DB_Y_Z = index_DB_Y_Z, quanti = convert.num,
                        nominal = nominal, ordinal = ordinal, logic = logic,
                        convert_num = convert.num, convert_class = convert.class,
                        prep_choice = "E")

  inst <- proxim_dist(dataB, norm = "E", prox = 0)
  nA <- inst$nA; nB <- inst$nB; Y <- inst$Y; Z <- inst$Z
  indXA <- inst$indXA; indXB <- inst$indXB; nbX <- length(indXA)
  Xobserv <- inst$Xobserv; Yobserv <- inst$Yobserv; Zobserv <- inst$Zobserv

  # --- 2. METRIC LEARNING (LDA) ---
  # Learning the optimal projection axis based on outcome classes in Database A
  message("Step 1: Learning Metric via LDA...")
  X_train <- as.data.frame(dataB[1:nA, -index_DB_Y_Z])
  Y_train <- as.factor(dataB[1:nA, index_DB_Y_Z[2]])
  lda_fit <- MASS::lda(X_train, Y_train)
  L <- lda_fit$scaling

  # Project all data into the learned space
  X_proj <- as.matrix(dataB[, -index_DB_Y_Z]) %*% L

  # --- 3. ROBUST CENTROID COST (M-ESTIMATOR) ---
  # Calculate centroids of outcome levels in the LDA space
  centroids_Y <- aggregate(X_proj[1:nA, ], list(dataB[1:nA, index_DB_Y_Z[2]]), mean)[,-1]
  centroids_Z <- aggregate(X_proj[(nA+1):(nA+nB), ], list(dataB[(nA+1):(nA+nB), index_DB_Y_Z[3]]), mean)[,-1]

  C_raw <- as.matrix(proxy::dist(centroids_Y, centroids_Z, method = "euclidean"))

  # Apply Huber M-estimator to dampen outliers in the cost matrix
  C_mat <- ifelse(abs(C_raw) <= k_huber,
                  0.5 * (C_raw^2),
                  k_huber * abs(C_raw) - 0.5 * k_huber^2)



  # --- 4. DUAL KL SINKHORN ENGINE ---
  # This solves the KL-unbalanced objective iteratively
  sinkhorn_dual_kl <- function(C, a, b, tA, tB, eps, iter = 100) {
    K <- exp(-C / eps)
    u <- rep(1, length(a)); v <- rep(1, length(b))
    fi_A <- tA / (tA + eps)
    fi_B <- tB / (tB + eps)
    for (i in 1:iter) {
      u <- (a / (as.vector(K %*% v) + 1e-15))^fi_A
      v <- (b / (t(K) %*% u + 1e-15))^fi_B
    }
    return(diag(as.vector(u)) %*% K %*% diag(as.vector(v)))
  }

  # --- 5. CORE TRANSPORT LOOP ---
  estimatorZA <- array(1/length(Z), dim = c(nbX, length(Y), length(Z)))
  estimatorYB <- array(1/length(Y), dim = c(nbX, length(Z), length(Y)))
  gamma_list <- list()

  for (x in 1:nbX) {
    mu_x <- sapply(Y, function(y) length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA) + 1e-10
    nu_x <- sapply(Z, function(z) length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB) + 1e-10

    gamma_x <- sinkhorn_dual_kl(C_mat, mu_x, nu_x, tau_A, tau_B, epsilon)
    gamma_list[[x]] <- gamma_x

    if (which.DB %in% c("A", "BOTH")) {
      for(y in 1:length(Y)) {
        rs <- sum(gamma_x[y, ])
        if(rs > 1e-12) estimatorZA[x, y, ] <- gamma_x[y, ] / rs
      }
    }
    if (which.DB %in% c("B", "BOTH")) {
      for(z in 1:length(Z)) {
        cs <- sum(gamma_x[, z])
        if(cs > 1e-12) estimatorYB[x, z, ] <- gamma_x[, z] / cs
      }
    }
  }

  # --- 6. INDIVIDUAL PREDICTIONS ---
  DATA1_OT <- dataB[dataB[, 1] == unique(dataB[, 1])[1], ]
  if (which.DB %in% c("A", "BOTH")) {
    probZA <- matrix(0, nA, length(Z))
    for(x in 1:nbX) for(i in indXA[[x]]) probZA[i, ] <- estimatorZA[x, Yobserv[i], ]
    DATA1_OT$OTpred <- factor(levels(dataB[, 3])[apply(probZA, 1, which.max)], levels = levels(dataB[, 3]))
  }

  DATA2_OT <- dataB[dataB[, 1] == unique(dataB[, 1])[2], ]
  if (which.DB %in% c("B", "BOTH")) {
    probYB <- matrix(0, nB, length(Y))
    for(x in 1:nbX) for(i in indXB[[x]]) probYB[i, ] <- estimatorYB[x, Zobserv[i + nA], ]
    DATA2_OT$OTpred <- factor(levels(dataB[, 2])[apply(probYB, 1, which.max)], levels = levels(dataB[, 2]))
  }

  tend <- Sys.time()

  # --- 7. OUTPUT ---
  res <- list(time_exe = difftime(tend, tstart),
              gamma = Reduce("+", gamma_list),
              profile = data.frame(ID = paste0("P_", 1:nrow(unique(Xobserv))), unique(Xobserv)),
              res_prox = inst,
              estimatorZA = estimatorZA,
              estimatorYB = estimatorYB,
              DATA1_OT = DATA1_OT,
              DATA2_OT = DATA2_OT,
              learned_metric = L)
  class(res) <- "otres"
  return(res)
}