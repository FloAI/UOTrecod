OT_LDA_unbalanced_joint <- function(datab, index_DB_Y_Z = 1:3,
                                    nominal = NULL, ordinal = NULL, logic = NULL,
                                    convert.num = NULL, convert.class = NULL,
                                    tau = 1.0, epsilon = 0.05,
                                    which.DB = "BOTH") {

  tstart <- Sys.time()

  # 1. Pre-processing
  dataB <- transfo_dist(datab, index_DB_Y_Z = index_DB_Y_Z, quanti = convert.num,
                        nominal = nominal, ordinal = ordinal, logic = logic,
                        convert_num = convert.num, convert_class = convert.class,
                        prep_choice = "E")

  inst <- proxim_dist(dataB, norm = "E", prox = 0)
  nA <- inst$nA; nB <- inst$nB; Y <- inst$Y; Z <- inst$Z
  indXA <- inst$indXA; indXB <- inst$indXB; nbX <- length(indXA)
  Xobserv <- inst$Xobserv; Yobserv <- inst$Yobserv; Zobserv <- inst$Zobserv

  # 2. Metric Learning via LDA
  message("Learning optimal projection axes via LDA...")
  X_train <- as.data.frame(dataB[1:nA, -index_DB_Y_Z])
  Y_train <- as.factor(dataB[1:nA, index_DB_Y_Z[2]])

  # Linear Discriminant Analysis
  lda_fit <- MASS::lda(X_train, Y_train)
  L <- lda_fit$scaling

  # 3. Calculate Centroid-based Cost Matrix using LDA Metric
  # We need the cost between levels of Y and levels of Z
  # We calculate the average X-profile (projected) for each level
  X_projected <- as.matrix(dataB[, -index_DB_Y_Z]) %*% L

  # Centroids for Y levels in Database A
  centroids_Y <- aggregate(X_projected[1:nA, ], list(dataB[1:nA, index_DB_Y_Z[2]]), mean)[,-1]
  # Centroids for Z levels in Database B
  centroids_Z <- aggregate(X_projected[(nA+1):(nA+nB), ], list(dataB[(nA+1):(nA+nB), index_DB_Y_Z[3]]), mean)[,-1]

  # C_mat is the distance between these outcome centroids in the LDA space
  C_mat <- as.matrix(proxy::dist(centroids_Y, centroids_Z, method = "euclidean"))

  # 4. Internal Sinkhorn Engine (KL-Unbalanced)
  sinkhorn_unbalanced <- function(C, a, b, tau, eps, iter = 60) {
    K <- exp(-C / eps)
    u <- rep(1, length(a)); v <- rep(1, length(b))
    fi <- tau / (tau + eps)
    for (i in 1:iter) {
      # u and v updates based on the KL divergence proximal operator
      u <- (a / (as.vector(K %*% v) + 1e-15))^fi
      v <- (b / (as.vector(t(K) %*% u) + 1e-15))^fi
    }
    return(diag(as.vector(u)) %*% K %*% diag(as.vector(v)))
  }

  # 5. Core Transport Loop
  estimatorZA <- array(1/length(Z), dim = c(nbX, length(Y), length(Z)))
  estimatorYB <- array(1/length(Y), dim = c(nbX, length(Z), length(Y)))
  gammaA_list <- list(); gammaB_list <- list()

  for (x in 1:nbX) {
    mu_x <- sapply(Y, function(y) length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA) + 1e-10
    nu_x <- sapply(Z, function(z) length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB) + 1e-10

    # mu_x must match rows of C_mat, nu_x must match columns
    gamma_x <- sinkhorn_unbalanced(C_mat, mu_x, nu_x, tau, epsilon)

    if (which.DB %in% c("A", "BOTH")) {
      for(y in 1:length(Y)) {
        row_s <- sum(gamma_x[y, ])
        if(row_s > 1e-12) estimatorZA[x, y, ] <- gamma_x[y, ] / row_s
      }
      gammaA_list[[x]] <- gamma_x
    }

    if (which.DB %in% c("B", "BOTH")) {
      for(z in 1:length(Z)) {
        col_s <- sum(gamma_x[, z])
        if(col_s > 1e-12) estimatorYB[x, z, ] <- gamma_x[, z] / col_s
      }
      gammaB_list[[x]] <- gamma_x
    }
  }

  # 6. Database Imputation
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

  res <- list(
    time_exe = difftime(tend, tstart),
    gamma_A = if(length(gammaA_list) > 0) Reduce("+", gammaA_list) else NULL,
    gamma_B = if(length(gammaB_list) > 0) Reduce("+", gammaB_list) else NULL,
    profile = data.frame(ID = paste0("P_", 1:nrow(unique(Xobserv))), unique(Xobserv)),
    res_prox = inst,
    estimatorZA = if(which.DB %in% c("A", "BOTH")) estimatorZA else NULL,
    estimatorYB = if(which.DB %in% c("B", "BOTH")) estimatorYB else NULL,
    DATA1_OT = DATA1_OT,
    DATA2_OT = DATA2_OT
  )
  class(res) <- "otres"
  return(res)
}