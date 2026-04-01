linear_our_log_posterior_reg_horseshoe <- function(Y, X, parameters, tau = 0.1, c = 1, sigma = 1) {
  # parameters: first element intercept, rest coefficients in order of columns of X
  beta0 <- parameters[1]
  beta <- parameters[-1]

  # prior for betas (horseshoe-marginal) and prior for intercept
  log_prior_vals <- sum(sapply(beta, function(bj) log_marginal_prior_reg_horseshoe(bj, tau, c)))
  log_prior_beta0 <- dnorm(beta0, mean = 0, sd = 10, log = TRUE)

  # Gaussian log-likelihood with known sigma
  n <- length(Y)
  mu <- beta0 + as.vector(X %*% beta)
  res <- Y - mu
  log_likelihood <- -0.5 * n * log(2 * pi * sigma^2) - sum(res^2) / (2 * sigma^2)

  return(log_prior_vals + log_prior_beta0 + log_likelihood)
}

linear_our_log_prior_reg_horseshoe <- function(Y, X, parameters, tau = 0.1, c = 1) {
  beta0 <- parameters[1]
  beta <- parameters[-1]
  log_prior_vals <- sum(sapply(beta, function(bj) log_marginal_prior_reg_horseshoe(bj, tau, c)))
  log_prior_beta0 <- dnorm(beta0, mean = 0, sd = 10, log = TRUE)
  return(log_prior_vals + log_prior_beta0)
}

linear_our_log_likelihood_reg_horseshoe <- function(Y, X, parameters, sigma = 1) {
  beta0 <- parameters[1]
  beta <- parameters[-1]
  n <- length(Y)
  mu <- beta0 + as.vector(X %*% beta)
  res <- Y - mu
  log_likelihood <- -0.5 * n * log(2 * pi * sigma^2) - sum(res^2) / (2 * sigma^2)
  return(log_likelihood)
}

# Main coordinate-descent function converted to linear regression
coordinate_descent_horseshoe_perm_XI_linear <- function(Y, X1, X2 = NULL, I_mat = NULL, Z = NULL,
                                                       tau = 0.1, c = 1, sigma = 1,
                                                       max_iter = 100, tol = 1e-3, K = 10,
                                                       verbose = TRUE) {

  # --- prepare matrices (same as your original function) --------------------
  if (is.null(X1)) stop("X1 (SNP matrix) must be provided.")
  X1 <- as.matrix(X1)
  n <- nrow(X1)
  p <- ncol(X1)                      # p = number of SNPs (NOT including X2)

  snp_names <- colnames(X1)
  if (is.null(snp_names)) snp_names <- paste0("snp", seq_len(p))

  if (!is.null(X2)) {
    Wmat <- as.matrix(X2)
    if (nrow(Wmat) != n) stop("X2 must have same number of rows as X1")
    if (is.null(colnames(Wmat))) colnames(Wmat) <- paste0("W", seq_len(ncol(Wmat)))
  } else {
    Wmat <- matrix(nrow = n, ncol = 0)
  }
  w <- ncol(Wmat)

  X_snp <- make_dom_rec_matrix(X1)

  if (!is.null(Z)) {
    Zmat <- as.matrix(Z)
    if (nrow(Zmat) != n) stop("Z must have same number of rows as X1")
    if (is.null(colnames(Zmat))) colnames(Zmat) <- paste0("cov", seq_len(ncol(Zmat)))
  } else {
    Zmat <- matrix(nrow = n, ncol = 0)
  }
  q <- ncol(Zmat)

  if (!is.null(I_mat)) {
    I_mat <- as.matrix(I_mat)
    if (nrow(I_mat) != n) stop("I_mat must have same number of rows as X1")
    r <- ncol(I_mat)
    if (is.null(colnames(I_mat))) colnames(I_mat) <- paste0("I", seq_len(r))
  } else {
    I_mat <- matrix(nrow = n, ncol = 0)
    r <- 0
  }

  X_new <- cbind(X_snp, Wmat, I_mat, Zmat)

  # --- helpers --------------------------------------------------------------
  softmax3 <- function(logv) {
    if (all(!is.finite(logv))) return(rep(NA_real_, length(logv)))
    m <- max(logv, na.rm = TRUE)
    exps <- exp(logv - m)
    exps / sum(exps)
  }

  safe_optimize <- function(obj, lower = -10, upper = 10, cur = 0, ngrid = 21) {
    try_grid <- seq(lower, upper, length.out = ngrid)
    vals_grid <- sapply(try_grid, function(x) {
      v <- tryCatch(obj(x), error = function(e) NA_real_, warning = function(w) NA_real_)
      if (!is.finite(v)) NA_real_ else v
    })
    if (all(is.na(vals_grid))) {
      return(list(success = FALSE, minimum = NA_real_, objective = NA_real_, reason = "no_finite_on_grid"))
    }
    opt_res <- tryCatch({
      withCallingHandlers({
        optimize(obj, c(lower, upper))
      }, warning = function(w) invokeRestart("muffleWarning"))
    }, error = function(e) NULL)
    if (!is.null(opt_res)) {
      val_at_opt <- tryCatch(obj(opt_res$minimum), error = function(e) NA_real_, warning = function(w) NA_real_)
      if (is.finite(val_at_opt)) {
        return(list(success = TRUE, minimum = opt_res$minimum, objective = val_at_opt, reason = "optimize_ok"))
      }
    }
    finite_idx <- which(is.finite(vals_grid))
    best_idx <- finite_idx[which.min(vals_grid[finite_idx])]
    best_x <- try_grid[best_idx]
    best_val <- vals_grid[best_idx]
    return(list(success = TRUE, minimum = best_x, objective = best_val, reason = "grid_fallback"))
  }

  # Gaussian log-likelihood helper (for incremental updates)
  log_lik <- function(Y, X_design, beta_full, sigma) {
    n <- length(Y)
    mu <- beta_full[1] + as.vector(X_design %*% beta_full[-1])
    res <- Y - mu
    - (0.5 * n * log(2 * pi * sigma^2) + sum(res^2) / (2 * sigma^2))
  }

  # This relies on linear_our_log_posterior_reg_horseshoe() defined above
  log_post_no_prior_on_Z <- function(Y, X_design, beta_full, tau, c, sigma) {
    if (q == 0) {
      return(linear_our_log_posterior_reg_horseshoe(Y, X_design, beta_full, tau, c, sigma))
    }
    beta_zeroZ <- beta_full
    z_start <- 1 + 3 * p + w + r + 1
    z_end <- 1 + 3 * p + w + r + q
    beta_zeroZ[z_start:z_end] <- 0

    lp_base <- linear_our_log_posterior_reg_horseshoe(Y, X_design, beta_zeroZ, tau, c, sigma)

    ll_full <- log_lik(Y, X_design, beta_full, sigma)
    ll_zeroZ <- log_lik(Y, X_design, beta_zeroZ, sigma)

    lp_base + (ll_full - ll_zeroZ)
  }
  # -----------------------------------------------------

  n_params <- 1 + ncol(X_new)
  best_beta <- rep(0, n_params)
  best_log_post <- -Inf
  final_perm_order <- NULL
  used_permutations <- character(0)

  best_snp_probs <- matrix(NA_real_, nrow = p, ncol = 3)
  colnames(best_snp_probs) <- c("ADD", "DOM", "REC")
  rownames(best_snp_probs) <- snp_names

  perm_count <- 0
  best_I_idx_perm_cols <- NULL

  while (perm_count < K) {
    permuted_snps <- sample.int(p)
    perm_str <- paste(permuted_snps, collapse = "-")
    if (perm_str %in% used_permutations) { next }
    used_permutations <- c(used_permutations, perm_str)
    perm_count <- perm_count + 1

    idx_perm_cols_snp <- unlist(lapply(permuted_snps,
                                       function(j) ((j - 1) * 3 + 1):((j - 1) * 3 + 3)))

    if (r > 0) {
      pref <- intersect(permuted_snps, seq_len(r))
      remaining <- setdiff(seq_len(r), pref)
      I_idx_perm_cols <- c(pref, remaining)
    } else {
      I_idx_perm_cols <- integer(0)
    }

    if (r > 0 && q > 0) {
      X_perm <- cbind(X_snp[, idx_perm_cols_snp, drop = FALSE],
                      Wmat,
                      I_mat[, I_idx_perm_cols, drop = FALSE],
                      Zmat)
    } else if (r > 0) {
      X_perm <- cbind(X_snp[, idx_perm_cols_snp, drop = FALSE],
                      Wmat,
                      I_mat[, I_idx_perm_cols, drop = FALSE])
    } else if (q > 0) {
      X_perm <- cbind(X_snp[, idx_perm_cols_snp, drop = FALSE],
                      Wmat,
                      Zmat)
    } else {
      X_perm <- cbind(X_snp[, idx_perm_cols_snp, drop = FALSE],
                      Wmat)
    }

    beta <- rep(0, 1 + ncol(X_perm))
    last_iter_mode_logps <- matrix(NA_real_, nrow = p, ncol = 3)

    for (iter in 1:max_iter) {
      beta_old <- beta

      # intercept
      obj_intercept <- function(b0) {
        tmp <- beta; tmp[1] <- b0
        -log_post_no_prior_on_Z(Y, X_perm, tmp, tau, c, sigma)
      }
      opt0 <- safe_optimize(obj_intercept, lower = -10, upper = 10, cur = beta[1])
      if (opt0$success) {
        beta[1] <- opt0$minimum
      } else {
        if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | intercept safe_optimize failed: ", opt0$reason)
      }

      # SNP blocks (3 columns each)
      for (block_pos in seq_len(p)) {
        block_cols <- ((block_pos - 1) * 3 + 1):((block_pos - 1) * 3 + 3)
        beta_indices <- 1 + block_cols

        real_name <- snp_names[permuted_snps[block_pos]]

        mode_logps <- rep(-Inf, 3)
        mode_vals <- rep(NA_real_, 3)
        mode_reasons <- rep(NA_character_, 3)

        for (mode in 1:3) {
          sel_idx <- beta_indices[mode]
          other_idxs <- beta_indices[-mode]

          obj_mode <- function(bsel) {
            tmp <- beta
            tmp[sel_idx] <- bsel
            tmp[other_idxs] <- 0
            -log_post_no_prior_on_Z(Y, X_perm, tmp, tau, c, sigma)
          }

          opt_mode <- safe_optimize(obj_mode, lower = -10, upper = 10, cur = beta[sel_idx])
          if (opt_mode$success) {
            tmp_beta <- beta
            tmp_beta[sel_idx] <- opt_mode$minimum
            tmp_beta[other_idxs] <- 0
            lp_mode <- log_post_no_prior_on_Z(Y, X_perm, tmp_beta, tau, c, sigma)
            mode_logps[mode] <- lp_mode
            mode_vals[mode] <- opt_mode$minimum
            mode_reasons[mode] <- opt_mode$reason

            if (opt_mode$reason == "grid_fallback" && verbose) {
              message("Perm ", perm_count, " | Iter ", iter, " | SNP ", real_name,
                      " | mode ", c("ADD","DOM","REC")[mode], " used grid fallback (block mode).")
            }
          } else {
            mode_logps[mode] <- -1e10
            mode_vals[mode] <- 0
            mode_reasons[mode] <- opt_mode$reason
            if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | SNP ", real_name,
                                 " | mode ", c("ADD","DOM","REC")[mode], " optimize failed. Reason: ", opt_mode$reason)
          }
        }

        last_iter_mode_logps[block_pos, ] <- mode_logps

        if (verbose) {
          probs_block <- softmax3(mode_logps)
          message("Perm ", perm_count, " | Iter ", iter,
                  " | SNP ", real_name,
                  " | ADD ", round(mode_logps[1], 4)," / P(ADD) = ", ifelse(is.na(probs_block[1]), "NA", sprintf("%.4f", probs_block[1])),
                  " | DOM ", round(mode_logps[2], 4)," / P(DOM) = ",  ifelse(is.na(probs_block[2]), "NA", sprintf("%.4f", probs_block[2])),
                  " | REC ", round(mode_logps[3], 4)," / P(REC) = ", ifelse(is.na(probs_block[3]), "NA", sprintf("%.4f", probs_block[3])),
                  " | choose = ", c("ADD","DOM","REC")[which.max(mode_logps)])
        }

        best_mode <- which.max(mode_logps)
        best_local_val <- mode_vals[best_mode]
        chosen_idx <- beta_indices[best_mode]
        other_chosen <- beta_indices[-best_mode]
        beta[chosen_idx] <- best_local_val
        beta[other_chosen] <- 0
      }

      # Update W (X2) coefficients
      if (w > 0) {
        W_start_in_Xperm <- 3 * p + 1
        W_end_in_Xperm <- 3 * p + w
        for (jj in W_start_in_Xperm:W_end_in_Xperm) {
          beta_idx <- 1 + jj
          obj_w <- function(bj) {
            tmp <- beta; tmp[beta_idx] <- bj
            -log_post_no_prior_on_Z(Y, X_perm, tmp, tau, c, sigma)
          }
          opt_w <- safe_optimize(obj_w, lower = -10, upper = 10, cur = beta[beta_idx])
          if (opt_w$success) {
            beta[beta_idx] <- opt_w$minimum
            if (opt_w$reason == "grid_fallback" && verbose) {
              message("Perm ", perm_count, " | Iter ", iter, " | updated W col at idx ", jj, " using grid fallback.")
            }
          } else {
            if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | W-block optimize failed for idx ", jj,
                                 " reason: ", opt_w$reason)
          }
        }
        if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | updated W block (", w, " cols).")
      }

      # Update I block
      if (r > 0) {
        I_start_in_Xperm <- 3 * p + w + 1
        I_end_in_Xperm <- I_start_in_Xperm + r - 1
        for (jj in I_start_in_Xperm:I_end_in_Xperm) {
          beta_idx <- 1 + jj
          obj_I <- function(bj) {
            tmp <- beta; tmp[beta_idx] <- bj
            -log_post_no_prior_on_Z(Y, X_perm, tmp, tau, c, sigma)
          }
          opt_i <- safe_optimize(obj_I, lower = -10, upper = 10, cur = beta[beta_idx])
          if (opt_i$success) {
            beta[beta_idx] <- opt_i$minimum
            if (opt_i$reason == "grid_fallback" && verbose) {
              message("Perm ", perm_count, " | Iter ", iter, " | updated I col at idx ", jj, " using grid fallback.")
            }
          } else {
            if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | I-block optimize failed for idx ", jj,
                                 " reason: ", opt_i$reason)
          }
        }
        if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | updated I block (", r, " cols).")
      }

      # Update Z covariates
      if (q > 0) {
        cov_start_in_Xperm <- 3 * p + w + r + 1
        cov_indices <- (1 + cov_start_in_Xperm):(1 + cov_start_in_Xperm + q - 1)
        for (j in cov_indices) {
          obj_cov <- function(bj) {
            tmp <- beta; tmp[j] <- bj
            -log_post_no_prior_on_Z(Y, X_perm, tmp, tau, c, sigma)
          }
          opt_c <- safe_optimize(obj_cov, lower = -10, upper = 10, cur = beta[j])
          if (opt_c$success) {
            beta[j] <- opt_c$minimum
            if (opt_c$reason == "grid_fallback" && verbose) {
              message("Perm ", perm_count, " | Iter ", iter, " | cov idx ", j, " used grid fallback (covariate).")
            }
          } else {
            if (verbose) message("Perm ", perm_count, " | Iter ", iter, " | cov optimize failed for idx ", j,
                                 " reason: ", opt_c$reason)
          }
        }
      }

      if (sqrt(sum((beta - beta_old)^2)) < tol) {
        if (verbose) message("Permutation ", perm_count, " converged at iteration ", iter)
        break
      }
    } # end iterations

    # Map back to full_beta in original (unpermuted) ordering:
    full_beta <- rep(NA_real_, 1 + 3 * p + w + r + q)
    full_beta[1] <- beta[1]

    for (orig in seq_len(p)) {
      pos <- which(permuted_snps == orig)
      perm_cols <- ((pos - 1) * 3 + 1):((pos - 1) * 3 + 3)
      orig_cols <- ((orig - 1) * 3 + 1):((orig - 1) * 3 + 3)
      full_beta[1 + orig_cols] <- beta[1 + perm_cols]
    }

    if (w > 0) {
      w_cols_perm_idx <- (3 * p + 1):(3 * p + w)
      full_beta[1 + w_cols_perm_idx] <- beta[1 + w_cols_perm_idx]
    }

    if (r > 0) {
      for (perm_pos in seq_len(r)) {
        orig_I_col <- I_idx_perm_cols[perm_pos]
        full_beta[1 + 3 * p + w + orig_I_col] <- beta[1 + 3 * p + w + perm_pos]
      }
    }

    if (q > 0) {
      from_idx <- (1 + 3 * p + w + r + 1)
      to_idx <- (1 + 3 * p + w + r + q)
      full_beta[from_idx:to_idx] <- beta[from_idx:to_idx]
    }

    # compute log posterior for this full_beta
    lp_now <- log_post_no_prior_on_Z(Y, X_new, full_beta, tau, c, sigma)

    # compute per-SNP model probs under full design
    snp_probs_now <- matrix(NA_real_, nrow = p, ncol = 3)
    colnames(snp_probs_now) <- c("ADD", "DOM", "REC")
    rownames(snp_probs_now) <- snp_names
    for (i in seq_len(p)) {
      orig_cols <- ((i - 1) * 3 + 1):((i - 1) * 3 + 3)
      logps <- rep(-Inf, 3)
      for (kcoord in 1:3) {
        idx_in_full <- 1 + orig_cols[kcoord]
        obj_k <- function(bk) {
          tmp <- full_beta
          tmp[idx_in_full] <- bk
          -log_post_no_prior_on_Z(Y, X_new, tmp, tau, c, sigma)
        }
        opt_k <- safe_optimize(obj_k, lower = -100, upper = 100, cur = full_beta[idx_in_full])
        if (opt_k$success) {
          tmp2 <- full_beta
          tmp2[idx_in_full] <- opt_k$minimum
          logps[kcoord] <- log_post_no_prior_on_Z(Y, X_new, tmp2, tau, c, sigma)
        } else {
          logps[kcoord] <- -1e10
        }
      }
      snp_probs_now[i, ] <- softmax3(logps)
    }

    last_mode_logps_orig_order <- matrix(NA_real_, nrow = p, ncol = 3)
    rownames(last_mode_logps_orig_order) <- snp_names
    colnames(last_mode_logps_orig_order) <- c("ADD","DOM","REC")
    for (block_pos in seq_len(p)) {
      orig_index <- permuted_snps[block_pos]
      last_mode_logps_orig_order[orig_index, ] <- last_iter_mode_logps[block_pos, ]
    }
    last_iter_probs_orig_order <- t(apply(last_mode_logps_orig_order, 1, softmax3))
    rownames(last_iter_probs_orig_order) <- snp_names
    colnames(last_iter_probs_orig_order) <- c("ADD","DOM","REC")

    if (verbose) {
      for (i in seq_len(p)) {
        logps_i <- last_mode_logps_orig_order[i, ]
        probs_i <- last_iter_probs_orig_order[i, ]
        message("Perm ", perm_count, " | SNP ", snp_names[i],
                " last-iter logps: ", paste0(sprintf("%.6g", logps_i), collapse = ", "),
                " -> probs (from last interaction): ", paste0(sprintf("%.4f", probs_i), collapse = ", "))
      }
    }

    if (lp_now > best_log_post) {
      best_log_post <- lp_now
      best_beta <- full_beta
      final_perm_order <- permuted_snps
      best_snp_probs <- last_iter_probs_orig_order

      if (r > 0) {
        best_I_idx_perm_cols <- I_idx_perm_cols
      } else {
        best_I_idx_perm_cols <- NULL
      }
    }

  } # end permutations

  message("Best permutation indices: ", paste(final_perm_order, collapse = "|"))
  message("Best log-posterior: ", best_log_post)

  if (r > 0) {
    if (!is.null(best_I_idx_perm_cols)) {
      if (!is.null(colnames(I_mat))) {
        reordered_names <- colnames(I_mat)[best_I_idx_perm_cols]
        message("I_mat column order in the BEST permutation (left-to-right in X_perm): ",
                paste0(reordered_names, collapse = "|"))
      } else {
        message("I_mat column order in the BEST permutation (left-to-right in X_perm) - original indices: ",
                paste0(best_I_idx_perm_cols, collapse = "|"))
      }
      pos_in_perm <- match(seq_len(r), best_I_idx_perm_cols)
      if (!is.null(colnames(I_mat))) {
        mapping_str <- paste0(colnames(I_mat), " -> pos_in_X_perm:", pos_in_perm, collapse = "; ")
        message("Mapping (original I_mat column -> position in X_perm): ", mapping_str)
      } else {
        mapping_str <- paste0(seq_len(r), " -> pos_in_X_perm:", pos_in_perm, collapse = "; ")
        message("Mapping (original I_mat index -> position in X_perm): ", mapping_str)
      }
    } else {
      message("I_mat was provided but no reordering information was saved for the best permutation.")
    }
  }

  names(best_beta) <- c("intercept",
                        colnames(X_snp),
                        if (w>0) colnames(Wmat) else character(0),
                        if (r>0) colnames(I_mat) else character(0),
                        if (q>0) colnames(Zmat) else character(0))

  message("Per-SNP posterior probabilities (ADD | DOM | REC):")
  for (i in seq_len(p)) {
    probs_i <- best_snp_probs[i, ]
    message(sprintf("%s : ADD = %.4f | DOM = %.4f | REC = %.4f",
                    snp_names[i],
                    probs_i[1], probs_i[2], probs_i[3]))
  }

  snp_probs_df <- as.data.frame(best_snp_probs, stringsAsFactors = FALSE)
  snp_probs_df$snp <- rownames(snp_probs_df)
  snp_probs_df <- snp_probs_df[, c("snp", "ADD", "DOM", "REC")]

  return(list(coefficients = best_beta,
              snp_model_probs = snp_probs_df,
              best_log_posterior = best_log_post,
              best_permutation = final_perm_order,
              best_I_idx_perm_cols = best_I_idx_perm_cols))
}
