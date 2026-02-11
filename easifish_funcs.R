# easifish functions

# here we load in a bunch of basic functions to keep the easifish code cleaner
# This includes some tidy functions, life history, fishery, selectivity functions etx
# Author Nick Hill

# Get productivity equation ----
# This functions estimates the productivity portion of EASIFish
# It  applies a length-based yield-per-recruit and spawn-per-recruit analysis.
# It returns a df by sim, spp, fishing mortality (x), and results including ypr, spr, slope and SSBr.
# The input data is a tible of format sim and spp column followed by a nested df.
# In the nested df is length, weight, maturity, delta-T, M and Fsusc.
get_Bprod3 <- function(nested_df, x_vals = seq(0, 10, 0.02)) {
  
  results <- vector("list", nrow(nested_df))
  
  for (i in seq_len(nrow(nested_df))) {
    
    sim <- nested_df$sim[i]
    spp <- nested_df$spp[i]
    df  <- nested_df$data[[i]][order(nested_df$data[[i]]$mn_length), ]
    
    nL <- nrow(df)
    nX <- length(x_vals)
    
    # ---- Pre-allocate matrices ----
    Fsusc_x <- outer(df$Fsusc, x_vals, "*")
    ZdT     <- (Fsusc_x + df$M) * df$dT
    
    surv <- exp(-apply(ZdT, 2, function(z)
      cumsum(c(0, z[-length(z)]))
    ))
    
    yprF <- df$wt *
      (Fsusc_x / (Fsusc_x + df$M)) *
      ((1 - exp(-ZdT)) * surv)
    
    SSBr <- df$wt * df$mat * surv
    
    # ---- Collapse across length bins ----
    ypr  <- colSums(yprF, na.rm = TRUE)
    SSBr <- colSums(SSBr, na.rm = TRUE)
    
    out <- data.frame(
      sim  = sim,
      spp  = spp,
      x    = x_vals,
      ypr  = ypr,
      SSBr = SSBr
    )
    
    out$slope <- c(diff(out$ypr) / diff(out$x), NA)
    out$spr   <- out$SSBr / out$SSBr[out$x == 0]
    
    results[[i]] <- out
  }
  
  dplyr::bind_rows(results)
}

# Get easi stats function ----
# This function extracts various results from EASI-Fish.
# It extracts fishing mortality and biomass related estimates and assoc. reference points.
# It takes an output of results from both the susceptibility and productivity portion of EASi-Fish.
get_easi_stats2 <- function(Bprod_summ, Fsusc_fshry) {
  # Compute F1 slope value
  # F1_value <- Bprod_summ |>
  #   dplyr::filter(x <= 0.02) |>
  #   group_by(sim, spp) |>
  #   mutate(slope = lag(((lead(ypr) - ypr) / (lead(x) - x)) * 0.1)) |>
  #   dplyr::filter(x == 0.02) |>
  #   dplyr::select(sim, spp, F1_slope = slope)
  # Get estimate of current fishing mortality
  Fest_vals <- Fsusc_fshry |>
    group_by(sim, spp) |>
    summarise(Fest = sum(Fage))
  
  get_F_at_spr_threshold <- function(df, threshold) {
    idx <- which.max(df$spr <= threshold)
    df$x[idx]
  }
  get_F_at_spr_threshold <- function(df, threshold, fallback = 5) {
    idxs <- which(df$spr >= threshold)
    if(length(idxs) == 0) {
      return(fallback)
    }
    idx <- max(idxs) + 1
    if(idx > nrow(df)) {
      return(fallback)
    }
    df$x[idx]
  }
  
  
  get_SSB_at_F <- function(df, F_val) {
    if (is.na(F_val)) return(NA_real_)
    if (F_val == 0) return(0)
    idx <- which(df$x <= F_val)
    if (length(idx) == 0) return(NA_real_)
    df$SSBr[idx[which.max(df$x[idx])]]
  }
  
  # Try one sim-spp at a time
  easi_stats <-
    Bprod_summ |>
    group_by(sim, spp) |>
    tidyr::nest() |>
    left_join(Fest_vals, by = c("sim", "spp")) |>
    mutate(FMSY = map_dbl(data, ~ .x$x[which.max(.x$ypr)]),
           F01 = map_dbl(data, ~ {
             valid_idxs <- which(.x$slope >= first(.x$slope) * 0.1)
             idx <- if(length(valid_idxs)) max(valid_idxs) + 1 else NA_integer_
             if (is.na(idx) || idx > nrow(.x)) idx <- nrow(.x)
             .x$x[idx]}),
           # F01 = map_dbl(data, ~ {idx <- which.min(abs(.x$slope - first(.x$slope) * 0.1))
           #   .x$x[idx]}),
           F20 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.2)),
           F40 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.4)),
           F60 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.6)),
           F80 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.8)),
           SSB    = map2_dbl(data, Fest,  ~ get_SSB_at_F(.x, .y)),
           SSBMSY = map2_dbl(data, FMSY,  ~ get_SSB_at_F(.x, .y)),
           SSB01  = map2_dbl(data, F01,   ~ get_SSB_at_F(.x, .y)),
           SSB20  = map2_dbl(data, F20,   ~ get_SSB_at_F(.x, .y)),
           SSB40  = map2_dbl(data, F40,   ~ get_SSB_at_F(.x, .y)),
           SSB60  = map2_dbl(data, F60,   ~ get_SSB_at_F(.x, .y)),
           SSB80  = map2_dbl(data, F80,   ~ get_SSB_at_F(.x, .y)),
           FFMSY = Fest/FMSY, SSBSSBMSY = SSB/SSBMSY,
           FF01 = Fest/F01, SSBSSB40 = SSB/SSB40) |>
    dplyr::select(-data) |>
    pivot_longer(-c(sim, spp), names_to = 'param', values_to = 'value')
  
  return(easi_stats)
  
}

# susceptibility equation ----
# This function calculates susceptibility
# Need to define each param needed as a column.
# olap type = 'annual' or 'month' which determines to use original or updated olap approach.
get_Ffinite <- function(olap_type = 'annual', olap, season, avail, encount, sel, avm, prm) {
  # Ensure olap_type is valid
  if (!olap_type %in% c("month", "annual")) {
    stop("Invalid olap_type. Must be either 'month' or 'annual'.")
  }
  
  # Compute Ffinite based on olap_type
  if (olap_type == "month") {
    Ffinite <- (olap * encount * sel) * (avm + ((1 - avm) * prm))
  } else {  # olap_type == "annual"
    Ffinite <- (olap * season * avail * encount * sel) * (avm + ((1 - avm) * prm))
  }
  
  return(Ffinite)
}

# get_Fadj <- function(Ffinite, Q, effort) {
#   Fadj <- (Ffinite * Q) * effort
#   return(Fadj)
# }
# get_Fage <- function(Fadj) {
#   Fage <- -log(1 - Fadj)
#   return(Fage)
# }

# prep fishery data function ----
# Takes a fishery df and converts it to a list with each element a unique species.
# Allows you to apply variability to different parameters by simulation. 
# Requires vary_param2 and easi_param2 functions (below) to work
# eg: param_settings = list(olap = list(dist = "lognormal", cv = 0.1)) as example
# Distributions available include: 'fixed', 'normal', 'normal01',
# 'uniform', 'lognormal', 'truncnorm'.
# run easi_param2 for details
prep_fishery_data <- function(fish_df, param_settings, verbose = TRUE) {
  if (verbose) message("🎣 Preparing fishery simulations...")
  
  spp_list <- unique(fish_df$spp)
  spp_data <- vector("list", length = length(spp_list))
  names(spp_data) <- spp_list
  
  for (spp in spp_list) {
    if (verbose) message("🟢 Processing species: ", spp)
    
    spp_df <- fish_df %>% filter(spp == !!spp)
    fishery_data <- vector("list", length = nrow(spp_df))
    
    for (i in seq_len(nrow(spp_df))) {
      row <- spp_df[i, ]
      nsim <- row$nsim
      
      # Repeat row nsim times
      sim_data <- row[rep(1, nsim), ] %>%
        rowid_to_column("sim")
      
      # Apply vary_param to specified parameters
      for (param in names(param_settings)) {
        if (param %in% names(sim_data)) {
          sim_data[[param]] <- vary_param2(sim_data[[param]], name = param, param_settings = param_settings)
        }
      }
      
      # Drop 'nsim' column
      sim_data <- sim_data %>% ungroup |> select(-nsim)
      
      fishery_data[[i]] <- sim_data
    }
    
    spp_data[[spp]] <- bind_rows(fishery_data)
  }
  
  if (verbose) message("✅ Fishery simulations prepared.")
  return(spp_data)
}

# encounterability function ----
# this calculates species encounterability based on the vertical distribution of spp vs fishery
# As per Griffiths et al. 2019.
# Need min and max depth of species, and of fishing gear to calculate proportional overlap from 0-1.
encount <- function(min_spp, max_spp, min_gr, max_gr) {
  if (max_spp <= min_gr | min_spp >= max_gr) { 
    enc <- 0
  } else if (max_spp <= max_gr & min_spp >= min_gr) { 
    enc <- 1
  } else { 
    enc <- (max(min(max_gr, max_spp) - max(min_gr, min_spp), 0)) / (max_spp - min_spp)
  }
  return(enc)
}

# easi_param function ----
# This function takes in lh and fishery params and applies a distribution to them (fixed, normal, or uniform)
# Each param has a row for: estimate, distribution, min, max, sd, and se from input csv.
easi_param <- function(x, fsry = NA, pmtr = NA) { 
  if(is.na(fsry) || length(fsry) == 0) {
    y <- x
  } else {
    # If fshry is a character vector with more than one element,
    # select only those columns of x that match the elements of fsry
    if(length(fsry) > 1) {
      y <- x[fsry]
    } else {
      # If fshry is a single character string, select the corresponding column of x
      y <- x[[fsry]]
    }
  }
  nms <- c(pmtr, paste0(pmtr,'_dist'), paste0(pmtr,'_max'), paste0(pmtr,'_min'), paste0(pmtr,'_sd'), paste0(pmtr,'_se'))
  y <- y[nms] |> map(type.convert, as.is = TRUE)
  
  y <- y[order(names(y))]
  p = y[[1]]; dist = y[[2]]; max = y[[3]]; min = y[[4]]; sd = y[[5]]; if(length(y) == 6){se = y[[6]]}
  if(is.na(p) | p == ""){return(NA)}
  if(dist == 'fixed') {z <- p}
  if(dist == 'normal') {z <- rnorm(1,p,sd)}
  if(dist == 'uniform') {z <- runif(1,min,max)}
  return(as.numeric(z))
}

# easi_tidy function ---- 
# This function takes a df with first column param, and then following column to be spp or fleets
# and puts them into a nested list of length spp/fleet and length(params) within each element
# Used in fleet_data script
easi_tidy <- function(df) {
  params <- df$param
  nms <- names(df)[-1]  # Exclude the 'param' column
  y <- list()
  
  map(nms, ~ {y[[.x]] <<- setNames(as.list(df[[.x]]), params)})
  
  return(y)
}

# length-weight function ----
# Calculate weight from length.
# Check grams/kg conversion. Shane /1000 but leads to decimal values???
lw <- function(a,len,b) {
  wt <- (a*len^b)#/1000
  return(wt)
}

# Delta-T2 fucntion ----
# Written as in Griffiths et al 2023 report equation and chen and gordon 1997
# Estimates the time taken to shift from one growth increment to the next. Used in YPR.
deltaT2 <- function(K, Linf, len, len_int) {
  # Apply rules
  len_adj <- case_when(
    len == 0      ~ NA_real_,          # handle 0 separately below
    len > Linf    ~ Linf - 0.1,
    TRUE          ~ len)
  
  len2 <- len_adj - len_int
  
  dT <- ifelse(len == 0,0,
               (1 / K) * log((Linf - len2) / (Linf - len2 - len_int)))
  
  return(dT)
}

# Maturity function ----
# Calculate maturity using L50, r and/or Lm 
# Setup so that r can't be negative
# Applies Griffiths et al. 2019 flowchart of estimating maturity based on available data

easi_mat2 <- function(len, L50 = NA, r = NA, Lm = NA) {
  if (!all(is.na(L50)) && !all(is.na(r))) {
    # Logistic curve
    return(1 / (1 + exp(-r * (len - L50))))
  }
  
  if (!all(is.na(L50)) && all(is.na(r))) {
    # Knife-edge at L50
    return(ifelse(len < L50, 0, 1))
  }
  
  if (all(is.na(L50)) && all(is.na(r)) && !all(is.na(Lm))) {
    # Knife-edge at Lm
    return(ifelse(len < Lm, 0, 1))
  }
  
  # Default: unknown params
  return(rep(NA_real_, length(len)))
}

# Calculate growth ----
# easi_growth converts age from length across various equations
# Need to define relevant input lh parameters.
# Available growth functions include:
# SVBGF, GVBGF, VBL0, Gompertz, Gompertz2, Richards, Schnute, Schnute2,
# Logistic, Logistic2, Logistic3, Robertson
easi_growth2 <- function(model, len, Linf, K, t0 = NA, L0 = NA, D = NA, m = NA,
                         age1 = NA, age2 = NA, len1 = NA, len2 = NA,
                         a = NA, b = NA, alpha = NA, g = NA) {
  
  n <- length(len)
  
  # Helper to recycle scalars to vector length
  r <- function(x) if (length(x) == 1) rep(x, n) else x
  
  model <- r(model); Linf <- r(Linf); K <- r(K); t0 <- r(t0); L0 <- r(L0)
  D <- r(D); m <- r(m); age1 <- r(age1); age2 <- r(age2)
  len1 <- r(len1); len2 <- r(len2); a <- r(a); b <- r(b); alpha <- r(alpha); g <- r(g)
  
  out <- numeric(n)
  
  for (i in seq_len(n)) {
    out[i] <- switch(model[i],
                     "SVBGF"     = log(1 - (len[i]/Linf[i])) / -K[i] + t0[i],
                     "GVBGF"     = log(1 - (len[i]/Linf[i])^D[i]) / -K[i] + t0[i],
                     "VBL0"      = 1 / K[i] * log((Linf[i] - L0[i]) / (Linf[i] - len[i])),
                     "Gompertz"  = -log(-log(len[i] / Linf[i]) / (1/K[i])) / K[i] + t0[i],
                     "Gompertz2" = -1 / K[i] * log(1 - log(len[i] / L0[i]) / log(Linf[i] / L0[i])), # SPL - Drew et al. 2015
                     "Richards"  = -log(((Linf[i]/len[i])^m[i] - 1) / t0[i]) / K[i],
                     "Schnute"   = age1[i] - (1/K[i]) * log((len[i] - Linf[i]) / (len1[i] - Linf[i])), # ALV - Teo et al. 2018
                     "Schnute2"  = age1[i] + log(1 - ((len[i]^b[i] - len1[i]^b[i]) / 
                                                        (len2[i]^b[i] - len1[i]^b[i])) *
                                                   (1 - exp(-a[i] * (age2[i] - age1[i])))) / (-a[i]), # WSH - Natanson and Skomal 2015
                     "Logistic"  = (1 / K[i]) * log((Linf[i] - alpha[i]) / (Linf[i] - len[i])), # BRO - Drew et al. 2016
                     "Logistic2" = (1 / g[i]) * log((len[i] * L0[i] - len[i] * Linf[i]) / 
                                                      (len[i] * L0[i] - Linf[i] * L0[i])), # FAL - Grant et al. 2018
                     "Logistic3" = (1 / K[i]) * log((len[i] * L0[i] - len[i] * Linf[i]) / 
                                                      (len[i] * L0[i] - Linf[i] * L0[i])), # BLR - Chin et al. 2013
                     "Robertson" = (t0[i] - log((Linf[i] / len[i]) - 1)) / K[i], #DUS - Joung et al. 2015
                     NA_real_
    )
  }
  
  return(out)
}

# Prepare life history data ----
# This function converts raw lh df into a list by species and varies params where defined.
# Requires vary_param2 and easi_param2 functions (below) to work
# param_settings = list(K = list(dist = "lognormal", cv = 0.1)) as example
# Distributions available include: 'fixed', 'normal', 'normal01',
# 'uniform', 'lognormal', 'truncnorm'.
# run easi_param2 for details
prep_life_history <- function(lh_df, param_settings, verbose = TRUE) {
  
  
  if (verbose) message("📥 Preparing life history simulations...")
  
  # Get unique length intervals
  unique_len_ints <- lh_df %>%
    pull(len_int) %>%
    unique() %>%
    na.omit() %>%
    sort()
  
  if (verbose) message("🔢 Found length intervals: ", paste(unique_len_ints, collapse = ", "))
  
  # Precompute length bins for each len_int
  bins_list <- map(unique_len_ints, function(li) {
    max_Linf <- max(lh_df$Linf[lh_df$len_int == li], na.rm = TRUE)
    seq(0, ceiling((max_Linf + li) / li) * li, by = li)
  })
  names(bins_list) <- as.character(unique_len_ints)
  
  if (verbose) message("📊 Length bins precomputed")
  
  spp_list <- lh_df$spp
  spp_data <- vector("list", length = length(spp_list))
  names(spp_data) <- spp_list
  
  for (i in seq_along(spp_list)) {
    spp_row <- lh_df[i, ]
    spp <- spp_row$spp
    if (verbose) message("🟢 Processing species: ", spp)
    
    nsim <- as.integer(spp_row$nsim)
    li <- spp_row$len_int
    len_bins <- bins_list[[as.character(li)]]
    
    # Expand rows by simulation
    sim_params <- spp_row[rep(1, nsim), ] %>%
      rowid_to_column(var = 'sim')
    
    # Apply variation to parameters
    for (param in names(param_settings)) {
      if (param %in% names(sim_params)) {
        sim_params[[param]] <- vary_param2(sim_params[[param]], name = param, param_settings = param_settings)
      }
    }
    
    # Ensure numeric and default values
    sim_params <- sim_params %>%
      mutate(
        Linf = as.numeric(Linf),
        L0 = as.numeric(L0 %||% 0),
        len_int = li
      )
    
    # Expand by length
    expanded_df <- sim_params %>%
      distinct(spp, sim, L0, Linf, .keep_all = TRUE) %>%
      crossing(length = len_bins) %>%
      #filter(length >= 0, length <= (Linf + li)) %>%
      mutate(mn_length = length - (li / 2)) |>
      #mutate(mn_length = (length + dplyr::lag(length, default = 0)) / 2) |>
      mutate(spp = spp, max_age = max_age) |>
      dplyr::filter(mn_length >= 0 & mn_length <= (Linf + (+ li / 2)))
    
    # Clean columns
    expanded_df <- expanded_df %>%
      select(-matches("_dist|_sd|_se|_min|_max")) %>%
      select(spp, sim, length, mn_length, max_age, everything())
    
    spp_data[[spp]] <- expanded_df
  }
  
  if (verbose) message("✅ All species prepared.")
  return(spp_data)
}

# Vary input parameters ----
# This function varies params based on input values and a distribution
vary_param2 <- function(x, name, param_settings) {
  setting <- param_settings[[name]]
  
  if (is.null(setting)) return(x)
  
  x <- as.numeric(x)
  
  vapply(x, function(val) {
    # Evaluate min/max if they are functions of val
    min_val <- if (is.function(setting$min)) setting$min(val) else setting$min
    max_val <- if (is.function(setting$max)) setting$max(val) else setting$max
    
    easi_param2(
      est  = val,
      dist = setting$dist %||% "fixed",
      min  = min_val,
      max  = max_val,
      sd   = setting$sd,
      se   = setting$se,
      cv   = setting$cv
    )
  }, numeric(1))
}

#Vary input parameters ----
#This function varies params based on input values and a distribution
vary_param <- function(x, name) {
  setting <- param_settings[[name]]
  if (is.null(setting)) return(x)  # No variation for this param

  vapply(x, function(val) {
    easi_param2(est = val,
                dist = setting$dist %||% "fixed",
                max = setting$max,
                min = setting$min,
                sd  = setting$sd,
                se  = setting$se,
                cv  = setting$cv)
  }, numeric(1))
}

easi_param2 <- function(est, dist = 'fixed', max = NA, min = NA, sd = NA, se = NA, cv = NA) {
  library(truncnorm)
  if (is.na(est) || est == "") return(NA)
  
  if (dist == 'fixed') {
    z <- est
    
  } else if (dist == 'normal') {
    z <- rnorm(1, mean = est, sd = sd)
    
  } else if (dist == 'normal01') {
    z <- rnorm_clipped(est, sd, min = 0, max = 1)
    
  } else if (dist == 'uniform') {
    z <- runif(1, min = min, max = max)
    
  } else if (dist == 'lognormal') {
    if (is.na(cv)) stop("Must provide cv for lognormal")
    sdlog <- sqrt(log(1 + cv^2))
    meanlog <- log(est) - 0.5 * sdlog^2
    z <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    
  } else if (dist == 'truncnorm') {
    if(is.na(sd)) stop("missing SD for truncnorm")
    z <- rtruncnorm(1, a=min, b=max, mean=est, sd=sd)
  } 
  else {
    z <- NA
  }
  
  return(z)
}

# Estimating natural mortality ----
# This function reads in predefined natural mortality functions (below) 
# and provides estimates based on lh input params
easiM_estimate <- function(funcs, ...) {
  args <- list(...)
  
  # Get individual M estimates
  M_vals <- sapply(funcs, function(f_name) {
    f <- tryCatch(match.fun(f_name), error = function(e) return(NA))
    if (is.function(f)) {
      f_args <- names(formals(f))
      usable_args <- args[names(args) %in% f_args]
      result <- tryCatch(do.call(f, usable_args), error = function(e) NA)
      return(as.numeric(result))  # Ensure the result is numeric
    } else {
      return(NA)
    }
  })
  
  mean_M <- as.vector(mean(M_vals, na.rm = TRUE))
  
  # Return as a data frame
  return(mean_M)
}

# Apply a range of natural mortaility functions and extract mean and 95% confidence bounds
easiM_estimator <- function(funcs, CV = 0.2, ...) {
  args <- list(...)
  
  # Get individual M estimates
  M_vals <- sapply(funcs, function(f_name) {
    f <- tryCatch(match.fun(f_name), error = function(e) return(NA))
    if (is.function(f)) {
      f_args <- names(formals(f))
      usable_args <- args[names(args) %in% f_args]
      tryCatch(do.call(f, usable_args), error = function(e) NA)
    } else {
      NA
    }
  })
  
  mean_M <- mean(M_vals, na.rm = TRUE)
  
  # Handle edge case
  if (is.na(mean_M) || mean_M <= 0) {
    return(data.frame(mean_M = NA, lower = NA, upper = NA))
  }
  
  # Lognormal distribution parameters from mean and CV
  sdlog <- sqrt(log(1 + CV^2))
  meanlog <- log(mean_M) - (sdlog^2) / 2
  
  # 95% confidence interval from lognormal
  lower <- qlnorm(0.025, meanlog, sdlog)
  upper <- qlnorm(0.975, meanlog, sdlog)
  
  # Return as a data frame
  return(data.frame(M_est = mean_M, M_lower = lower, M_upper = upper))
}

# easiM_estimator(funcs = c("hoenig83", "Liu_tmax", "Liu_Linf",
#                            "jensen", "pauly_sst","pauly_kt", "durueil",
#                            "durueil_shark", "hoenig_fish", "durueil_Linf", "frisk1",
#                            "Then_tmax", "Then_Linf"), CV = 0.1, tmax = 20, K = 0.1, L0 = 60)

# easiM_estimator_all function ----
# This function reads in predefined natural mortality functions (below) 
# and provides estimates based on lh input params
# Same function as above, but returns a df of all M values rather than a mean with bounds, plus optional plot

easiM_estimator_all <- function(funcs, plot = FALSE, ...) {
  args <- list(...)
  
  results <- lapply(funcs, function(f_name) {
    f <- tryCatch(match.fun(f_name), error = function(e) return(NULL))
    if (!is.null(f)) {
      f_args <- names(formals(f))
      usable_args <- args[names(args) %in% f_args]
      
      M_val <- tryCatch(do.call(f, usable_args), error = function(e) NA)
      
      if (!is.na(M_val)) {
        return(data.frame(func_name = f_name, M = M_val))
      }
    }
    return(NULL)
  })
  
  # Combine valid rows into a data frame
  M_df <- do.call(rbind, results)
  
  # Add row with mean M
  if (!is.null(M_df) && nrow(M_df) > 0) {
    mean_M <- mean(M_df$M, na.rm = TRUE)
    M_df <- rbind(M_df, data.frame(func_name = "mean", M = mean_M))
  } else {
    M_df <- data.frame(Estimator = "mean", M_est = NA)
  }
  # Optional plot
  if (plot) {
    p <- ggplot(M_df[M_df$func_name != "mean", ], aes(x = func_name, y = M)) +
      geom_point(color = "blue", size = 3) +
      geom_hline(yintercept = M_df$M[M_df$func_name == "mean"], linetype = "dashed", color = "red") +
      labs(x = "Estimator", y = "Natural mortality") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0,NA)
    
    print(p)
  }
  return(M_df)
  
}

# easiM_estimator_all(funcs = c("hoenig83", "Liu_tmax", "Liu_Linf",
#                            "jensen", "pauly_sst","pauly_kt", "durueil",
#                            "durueil_shark", "hoenig_fish", "durueil_Linf", "frisk1",
#                            "Then_tmax", "Then_Linf"), CV = 0.1, tmax = 20, K = 0.1, L0 = 60)

# Lmax to Linf function ----
# Convert an Lmax estimate to Linf from Froese and Binolahn 2000
Lmax_to_Linf <- function(Lmax) {
  Linf <- 10 ^ (0.044 + 0.9841 * log(Lmax, base = 10))
  return(Linf)
}

# Hoenig 1983 M function ----
# Estimate M function from Hoenig 1983
hoenig83 <- function(tmax) {
  M <- 4.3/tmax
  return(M)
}

# Hoenig M function for fishes ----
hoenig_fish <- function(tmax) {
  M <- exp(1.46 - 1.01 * log(tmax))
  return(M)
}

# Pauly M function with Linf and SST ----
pauly_lt <- function(Linf, sst, Lmax) {
  if(!is.na(Linf)) {M <- 10 ^ (0.556 - 0.718 * log(Linf)) + (0.02 * sst)}
  if(is.na(Linf)) {Linf <- 10 ^ (0.044 + 0.9841 * log(Lmax, base = 10))
  M <- 10 ^ (0.556 - 0.718 * log(Linf)) + (0.02 * sst)}
  return(M)
}

# Pauly M function with Linf, K and sst ----
pauly_sst <- function(sst, Linf, K) {
  if (sst < 0.001 || is.na(Linf) || is.na(K) || Linf == "" || K == "") {
    return("")
  } else {
    M <- exp(-0.0066 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(sst))
    return(M)
  }
}

# Pauly M function with K and sst ----
pauly_kt <- function(sst, K) {
  if (is.na(sst) || is.na(K) || sst == "" || K == "") {
    return("")
  } else {
    M <- K * exp(-0.22 + 0.3 * log(sst))
    return(M)
  }
}

# jensen M function ----
jensen <- function(K) {
  if(is.na(K)) { return("")}
  else{
    M <- 1.6 * K
    return(M)
  }
}

# Durueil et al. 2021 M function for all species ----
durueil <- function(tmax) {
  if(is.na(tmax) | !is.numeric(tmax)) { return("Error: Check tmax")}
  M <- exp(1.551 - 1.061 * log(tmax))
  return(M)
}

# Durueil et al. 2021 M function for sharks ----
durueil_shark <- function(tmax) {
  if(is.na(tmax) | !is.numeric(tmax)) { return("Error: Check tmax")}
  M <- exp(1.583 - 1.087 * log(tmax))
  return(M)
}

# Durueil et al. 2021 M function with Linf ----
durueil_Linf <- function(Linf, L0, K, X = 0.99) {
  get_Tmax <- function(Linf, K, L0, X) {
    tmax <- (1/K) * log((Linf - L0)/((1 - X) * L0))
    return(tmax)
  }
  M <- exp(1.583 - 1.087 * log(get_Tmax(Linf, K, L0, X)))
  return(M)
}

# Then et al. 2015 M function with T max also known as Hoenig_nls ----
then_Tmax <- function(tmax) {
  M <- 4.899 * tmax ^ -0.916
  return(M)
}

# Then et al. 2015 M function with Linf also known as Pauly_nls ----
then_Linf <- function(K, Linf) {
  M <- 4.118 * (K ^ 0.73) * (Linf ^ -0.33)
  return(M)
}

# Frisk et al 2001 M function for sharks with K ----
frisk1 <- function(K) {
  M <- 0.42 * log(K) - 0.83
  M <- exp(M)
  return(M)
}

# Frisk et al 2001 M function for sharks with age at maturity ----
frisk2 <- function(tmat) {
  M <- 1 / (0.44 * tmat + 1.87)
  return(M)
}

# Liu_tmax function ----
# Liu et al 2021 natural mortality function using BRT
# Note, you need metadata R data file.
# get from here - https://github.com/ChanjuanLiu92/M-estimate

Liu_tmax <- function(K = NA, Linf = NA, tmax = NA, Liu_method = 'tmax'){
  library(caret)
  library(gbm)
  
  if(Liu_method == 'K') {Liu_method <- 'BRT1'}
  if(Liu_method == 'Linf') {Liu_method <- 'BRT2'}
  if(Liu_method == 'tmax') {Liu_method <- 'BRT3'}
  
  # Deifne data to predict to
  mydata <- data.frame(K = K, Linf = Linf, tmax = tmax, Class = 1)
  
  # Create data to model to
  load('metadata324.Rdata')
  M=metadata324$M
  K=metadata324$K
  Linf=metadata324$Linf
  tmax=metadata324$tmax
  Class=as.factor(metadata324$Class)
  
  data1=na.omit(data.frame(M,tmax,Class))
  data2=na.omit(data.frame(M,K,Linf,Class))
  data3=na.omit(data.frame(M,K,Linf,tmax,Class))
  
  if (Liu_method=="BRT1"){
    model1=gbm(M~.,data=data1,distribution="gaussian",
               n.trees=400,
               interaction.depth=5,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model1,mydata,n.trees=400))[1,]
  }
  
  else if (Liu_method=="BRT2"){
    model2=gbm(M~.,data=data2,distribution="gaussian",
               n.trees=300,
               interaction.depth=8,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model2,mydata,n.trees=300))[1,]
  }
  
  else if (Liu_method=="BRT3"){
    model3=gbm(M~.,data=data3,distribution="gaussian",
               n.trees=380,
               interaction.depth=8,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model3,mydata,n.trees=380))[1,]
  }
  return(out)
}

# Liu_Linf function ----
# Liu et al 2021 natural mortality function using BRT
# Note, you need metadata R data file.
# get from here - https://github.com/ChanjuanLiu92/M-estimate

Liu_Linf <- function(K = NA, Linf = NA, tmax = NA, Liu_method = 'Linf'){
  library(caret)
  library(gbm)
  
  if(Liu_method == 'K') {Liu_method <- 'BRT1'}
  if(Liu_method == 'Linf') {Liu_method <- 'BRT2'}
  if(Liu_method == 'tmax') {Liu_method <- 'BRT3'}
  
  # Define data to predict to
  mydata <- data.frame(K = K, Linf = Linf, tmax = tmax, Class = 1)
  
  # Create data to model to
  load('metadata324.Rdata')
  M=metadata324$M
  K=metadata324$K
  Linf=metadata324$Linf
  tmax=metadata324$tmax
  Class=as.factor(metadata324$Class)
  
  data1=na.omit(data.frame(M,tmax,Class))
  data2=na.omit(data.frame(M,K,Linf,Class))
  data3=na.omit(data.frame(M,K,Linf,tmax,Class))
  
  if (Liu_method=="BRT1"){
    model1=gbm(M~.,data=data1,distribution="gaussian",
               n.trees=400,
               interaction.depth=5,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model1,mydata,n.trees=400))[1,]
  }
  
  else if (Liu_method=="BRT2"){
    model2=gbm(M~.,data=data2,distribution="gaussian",
               n.trees=300,
               interaction.depth=8,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model2,mydata,n.trees=300))[1,]
  }
  
  else if (Liu_method=="BRT3"){
    model3=gbm(M~.,data=data3,distribution="gaussian",
               n.trees=380,
               interaction.depth=8,
               shrinkage=0.01,
               n.minobsinnode = 1)
    out <- data.frame(predict(model3,mydata,n.trees=380))[1,]
  }
  return(out)
}


























#####################################################
# Generic functions ----

# get_Bprod function ----
# This function runs the productivity component of the EASIFish assessment
# This includes ypr and spr
# Takes a nested df by spp and sim with a nested data df that has
# mn_length, wt, mat, dT, M, Fsusc, Fsusc_min

# get_Bprod2 <- function(nested_df, x_vals = seq(0, 5, 0.02)) {
#   results <- list()
#   
#   for (i in seq_len(nrow(nested_df))) {
#     sim <- nested_df$sim[i]
#     spp <- nested_df$spp[i]
#     df <- nested_df$data[[i]]
#     
#     df <- df[order(df$mn_length), ]
#     Fsusc_min_val <- unique(df$Fsusc_min)
#     out <- data.frame()
#     
#     
#     
#     for (x in x_vals) {
#       # Total mortality of length bin (whole population) - Z = M when F = 0.
#       ZdT <- ((df$Fsusc * x) + df$M) * df$dT
#       
#       # Total survivorship up to length bin
#       surv <- exp(-cumsum(dplyr::lag(ZdT, n=1, default = 0)))
#       
#       # YPR for fished portion of population
#       yprF <- (df$wt * (df$Fsusc * x) / ((df$Fsusc * x) + df$M)) *
#                        ((1 - exp(-ZdT)) * surv)
#       
#       # SSBr is weight of length bin * maturity * how many have survived to contribute
#       SSBr <- df$wt * df$mat * surv
#       
#       #surv2 <- cumprod(surv)
#       
#       tmp <- data.frame(df, ZdT, yprF, surv, SSBr) |>
#         mutate(across(where(is.numeric), ~round(.,4)))
#       
#       # To get total ypr and spr sum across all length bins.
#       out <- rbind(out, data.frame(sim = sim, spp = spp, x = x,
#                                    #ypr = mean(yprF, na.rm = TRUE),
#                                    ypr = sum(yprF, na.rm = TRUE),
#                                    SSBr = sum(SSBr, na.rm = TRUE)))
#     }
#     
#     out$slope <- c(diff(out$ypr) / diff(out$x), NA)
#     out$spr <- out$SSBr / out$SSBr[out$x == 0]
#     
#     results[[i]] <- out
#   }
#   
#   dplyr::bind_rows(results)
# }
# 
# get_Bprod3 <- function(nested_df, x_vals = seq(0, 10, 0.02)) {
#   
#   results <- vector("list", nrow(nested_df))
#   
#   for (i in seq_len(nrow(nested_df))) {
#     
#     sim <- nested_df$sim[i]
#     spp <- nested_df$spp[i]
#     df  <- nested_df$data[[i]][order(nested_df$data[[i]]$mn_length), ]
#     
#     nL <- nrow(df)
#     nX <- length(x_vals)
#     
#     # ---- Pre-allocate matrices ----
#     Fsusc_x <- outer(df$Fsusc, x_vals, "*")
#     ZdT     <- (Fsusc_x + df$M) * df$dT
#     
#     surv <- exp(-apply(ZdT, 2, function(z)
#       cumsum(c(0, z[-length(z)]))
#     ))
#     
#     yprF <- df$wt *
#       (Fsusc_x / (Fsusc_x + df$M)) *
#       ((1 - exp(-ZdT)) * surv)
#     
#     SSBr <- df$wt * df$mat * surv
#     
#     # ---- Collapse across length bins ----
#     ypr  <- colSums(yprF, na.rm = TRUE)
#     SSBr <- colSums(SSBr, na.rm = TRUE)
#     
#     out <- data.frame(
#       sim  = sim,
#       spp  = spp,
#       x    = x_vals,
#       ypr  = ypr,
#       SSBr = SSBr
#     )
#     
#     out$slope <- c(diff(out$ypr) / diff(out$x), NA)
#     out$spr   <- out$SSBr / out$SSBr[out$x == 0]
#     
#     results[[i]] <- out
#   }
#   
#   dplyr::bind_rows(results)
# }
# 
# 
# #tst <- get_Bprod3(Bprod_ins, x_vals = seq(0, 10, 0.02))
# #tst2 <- get_Bprod2(Bprod_ins, x_vals = seq(0, 2, 0.1))
# 
# get_Bprod <- function(nested_df, x_vals = seq(0, 10, 0.02)) {
#   results <- list()
#   
#   for (i in seq_len(nrow(nested_df))) {
#     sim <- nested_df$sim[i]
#     spp <- nested_df$spp[i]
#     df <- nested_df$data[[i]]
#     
#     df <- df[order(df$mn_length), ]
#     Fsusc_min_val <- unique(df$Fsusc_min)
#     out <- data.frame()
#     
# 
#     
#     for (x in x_vals) {
#       ypr_termj <- ifelse(df$Fsusc == 0, 0, ((df$Fsusc * x) + df$M) * df$dT)
#       ypr_termj2 <- ((df$Fsusc * x) + df$M) * df$dT
#       ypr_termk <- dplyr::lag(ypr_termj, n=1, default = 0)
#       ypr_termk2 <- dplyr::lag(ypr_termj2, n=1, default = 0)
#       #surv <- exp(-cumsum(replace(ypr_term, is.na(ypr_term), 0)))
#       yprF <- ifelse(ypr_termj == 0, 0, (df$wt * (df$Fsusc * x) / ((df$Fsusc * x) + df$M)) * 
#                          ((1 - exp(-ypr_termj)) * exp(-cumsum(ypr_termk))))
#       yprF2 <- ifelse(df$Fsusc == 0, 0, (df$wt * (df$Fsusc * x) / ((df$Fsusc * x) + df$M)) * 
#                        ((1 - exp(-ypr_termj)) * exp(-cumsum(ypr_termk))))
#       
#       tmp <- data.frame(df, ypr_termj, ypr_termj2, ypr_termk, ypr_termk2, 
#                        SSBlc,SSBlc2, SSBr, SSBr2)
# 
#       
#       # df$ypr_term <- ifelse(df$Fsusc == 0, 0, ((df$Fsusc * x) + df$M) * df$dT)
#       # #surv <- exp(-cumsum(replace(ypr_term, is.na(ypr_term), 0)))
#       # df$yprF <- ifelse(
#       #   ypr_term == 0,
#       #   0,
#       #   (df$wt * (df$Fsusc * x)) / ((df$Fsusc * x) + df$M) *
#       #     (1 - exp(-ypr_term)) * exp(-cumsum(ypr_term))
#       # )
#       
#       # df$SSBlc <- ifelse(df$Fsusc == 0, 0,
#       #                 ifelse(df$mn_length == Fsusc_min_val, 1,
#       #                        ifelse(df$Fsusc > 0 & df$mn_length > Fsusc_min_val,
#       #                               exp(-c(NA, head(ypr_term, -1))),
#       #                               NA_real_)))
#       
#       #SSBlc = ifelse(df$Fsusc == 0, 0, exp(-ypr_termk))
#       SSBlc <- exp(-ypr_termk)
#       SSBlc2 <- exp(-ypr_termk2)
#       
#       
#       # df$SSBlc <- ifelse(df$Fsusc == 0, 0,
#       #                     exp(-c(NA, head(ypr_term, -1))))
#       
#       
#       # SSBr <- ifelse(df$Fsusc == 0, NA, df$wt * df$mat * cumprod2(SSBlc))
#       #SSBr <- df$wt * df$mat * cumprod2(SSBlc)
#       SSBr <- round(df$wt * df$mat * cumprod(SSBlc), 5)
#       SSBr2 <- round(df$wt * df$mat * cumprod(SSBlc2), 5)
#       
# 
#       
#       out <- rbind(out, data.frame(sim = sim, spp = spp, x = x,
#                                    #ypr = mean(yprF, na.rm = TRUE),
#                                    ypr = sum(yprF, na.rm = TRUE),
#                                    SSBr = sum(SSBr, na.rm = TRUE)))
#       
#       # z <- data.frame(df, yprj = ypr_termj, yprk = ypr_termk, yprF = yprF,
#       #                 SSBlc = SSBlc, SSBr = SSBr) |>
#       #   mutate(across(where(is.numeric), ~ round(., 9)))
#     }
#     
#     out$slope <- c(diff(out$ypr) / diff(out$x), NA)
#     out$spr <- out$SSBr / out$SSBr[out$x == 0]
#     
#     results[[i]] <- out
#   }
#   
#   dplyr::bind_rows(results)
# }
# 
# # out |> 
# #   ggplot() +
# #   aes(x, spr, col = sim) +
# #   geom_line() +
# #   xlim(0,2)
# 
# # get Ffinite function----
# # get_Ffinite <- function(olap, season, avail, encount, sel, avm, prm) {
# #   Ffinite <- (olap * season * avail * encount * sel) * (avm + ((1 - avm) * prm))
# #   return(Ffinite)
# # }
# 
# get_Ffinite <- function(olap_type = 'annual', olap, season, avail, encount, sel, avm, prm) {
#   # Ensure olap_type is valid
#   if (!olap_type %in% c("month", "annual")) {
#     stop("Invalid olap_type. Must be either 'month' or 'annual'.")
#   }
#   
#   # Compute Ffinite based on olap_type
#   if (olap_type == "month") {
#     Ffinite <- (olap * encount * sel) * (avm + ((1 - avm) * prm))
#   } else {  # olap_type == "annual"
#     Ffinite <- (olap * season * avail * encount * sel) * (avm + ((1 - avm) * prm))
#   }
#   
#   return(Ffinite)
# }
# get_Fadj <- function(Ffinite, Q, effort) {
#   Fadj <- (Ffinite * Q) * effort
#   return(Fadj)
# }
# get_Fage <- function(Fadj) {
#   Fage <- -log(1 - Fadj)
#   return(Fage)
# }
# 
# get_susceptibility <- function(df, nsim = 1, param_settings = list()) {
# 
#   
#   # Store results from all simulations
#   out <- vector("list", nsim)
#   
#   for (i in seq_len(nsim)) {
#     temp <- df
#     
#     # Apply variability via easi_param2 for specified parameters
#     for (param in names(param_settings)) {
#       if (param %in% names(temp)) {
#         temp[[param]] <- vary_param2(temp[[param]], param)
#       }
#     }
#     
#     # Calculate Ffinite, Fadj, Fage
#     Ffinite <- get_Ffinite(
#       olap    = temp$olap,
#       season  = temp$season,
#       avail   = temp$avail,
#       encount = temp$encount,
#       sel     = temp$sel,
#       avm     = temp$avm,
#       prm     = temp$prm
#     )
#     
#     Fadj <- get_Fadj(Ffinite, temp$Q, temp$effort)
#     Fage <- get_Fage(Fadj)
#     
#     # Return results
#     out[[i]] <- data.frame(
#       sim = i,
#       spp = temp$spp,
#       fishery = temp$fishery,
#       mn_length = temp$mn_length,
#       Ffinite = Ffinite,
#       Fadj = Fadj,
#       Fage = Fage)
#   }
#   
#   # Combine and return
#   do.call(rbind, out)
# }
# 
# # get_easi_stats_safe <- function(Bprod_summ, Fsusc_fshry) {
# #   # Compute F1 slope value
# #   F1_value <- Bprod_summ |>
# #     dplyr::filter(x <= 0.02) |>
# #     group_by(sim, spp) |>
# #     mutate(slope = lag(((lead(ypr) - ypr) / (lead(x) - x)) * 0.1)) |>
# #     dplyr::filter(x == 0.02) |>
# #     dplyr::select(sim, spp, F1_slope = slope)
# #   # Get estimate of current fishing mortality
# #   totF <- Fsusc_fshry |>
# #     group_by(sim, spp) |>
# #     summarise(Fest = sum(Fage), susc = 1 - exp(-sum(Fage)), .groups = "drop")
# #   # Try one sim-spp at a time
# #   safe_results <- Bprod_summ |>
# #     group_by(sim, spp) |>
# #     tidyr::nest() |>
# #     dplyr::left_join(F1_value, by = c("sim", "spp")) |>
# #     dplyr::left_join(totF, by = c("sim", "spp")) |>
# #     mutate(stats = purrr::pmap(
# #       list(data, F1_slope, Fest, FMSY = NA_real_, F01 = NA_real_),
# #       function(df, F1_slope, Fest, FMSY, F01) {
# #         tryCatch({
# #           FMSY <- df$x[which.max(df$ypr)]
# #           F01 <- df$x[which(df$slope <= F1_slope)[1]]
# #           Fsprs <- purrr::map_dbl(c(0.2, 0.4, 0.6, 0.8), function(thresh) {
# #             idx <- which.max(df$spr <= thresh)
# #             if (length(idx) == 0 || is.na(idx)) return(NA_real_)
# #             df$x[idx]
# #           })
# #           names(Fsprs) <- c("F20", "F40", "F60", "F80")
# #           get_SSB_at_F <- function(df, F_val) {
# #             if (is.na(F_val)) return(NA_real_)
# #             if (F_val == 0) return(0)
# #             idx <- which(df$x <= F_val)
# #             if (length(idx) == 0) return(NA_real_)
# #             df$SSBr[idx[which.max(df$x[idx])]]
# #           }
# #           SSB_vals <- purrr::map_dbl(
# #             c(Fest, FMSY, F01, Fsprs),
# #             ~ get_SSB_at_F(df, .x)
# #           )
# #           names(SSB_vals) <- c("SSB", "SSBMSY", "SSB01", "SSB20", "SSB40", "SSB60", "SSB80")
# #           FFMSY <- Fest / FMSY
# #           SSBSSBMSY <- SSB_vals["SSB"] / SSB_vals["SSBMSY"]
# #           tibble::tibble(
# #             param = c("FMSY", "F01", names(Fsprs), names(SSB_vals), "Fest", "FFMSY", "SSBSSBMSY"),
# #             value = c(FMSY, F01, Fsprs, SSB_vals, Fest, FFMSY, SSBSSBMSY)
# #           )
# #         }, error = function(e) {
# #           message("Skipping sim = ", unique(df$sim), ", spp = ", unique(df$spp), ": ", e$message)
# #           return(NULL)
# #         })
# #       }
# #     )) |>
# #     tidyr::unnest(cols = stats) |>
# #     dplyr::select(sim, spp, param, value)
# #   return(safe_results)
# # }
# 
# get_easi_stats2 <- function(Bprod_summ, Fsusc_fshry) {
#   # Compute F1 slope value
#   # F1_value <- Bprod_summ |>
#   #   dplyr::filter(x <= 0.02) |>
#   #   group_by(sim, spp) |>
#   #   mutate(slope = lag(((lead(ypr) - ypr) / (lead(x) - x)) * 0.1)) |>
#   #   dplyr::filter(x == 0.02) |>
#   #   dplyr::select(sim, spp, F1_slope = slope)
#   # Get estimate of current fishing mortality
#   Fest_vals <- Fsusc_fshry |>
#     group_by(sim, spp) |>
#     summarise(Fest = sum(Fage))
# 
#   get_F_at_spr_threshold <- function(df, threshold) {
#     idx <- which.max(df$spr <= threshold)
#     df$x[idx]
#   }
#   get_F_at_spr_threshold <- function(df, threshold, fallback = 5) {
#     idxs <- which(df$spr >= threshold)
#     if(length(idxs) == 0) {
#       return(fallback)
#     }
#     idx <- max(idxs) + 1
#     if(idx > nrow(df)) {
#       return(fallback)
#     }
#     df$x[idx]
#   }
# 
# 
#   get_SSB_at_F <- function(df, F_val) {
#     if (is.na(F_val)) return(NA_real_)
#     if (F_val == 0) return(0)
#     idx <- which(df$x <= F_val)
#     if (length(idx) == 0) return(NA_real_)
#     df$SSBr[idx[which.max(df$x[idx])]]
#   }
# 
#   # Try one sim-spp at a time
#   easi_stats <-
#     Bprod_summ |>
#     group_by(sim, spp) |>
#     tidyr::nest() |>
#     left_join(Fest_vals, by = c("sim", "spp")) |>
#     mutate(FMSY = map_dbl(data, ~ .x$x[which.max(.x$ypr)]),
#            F01 = map_dbl(data, ~ {
#              valid_idxs <- which(.x$slope >= first(.x$slope) * 0.1)
#              idx <- if(length(valid_idxs)) max(valid_idxs) + 1 else NA_integer_
#              if (is.na(idx) || idx > nrow(.x)) idx <- nrow(.x)
#              .x$x[idx]}),
#       # F01 = map_dbl(data, ~ {idx <- which.min(abs(.x$slope - first(.x$slope) * 0.1))
#       #   .x$x[idx]}),
#       F20 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.2)),
#       F40 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.4)),
#       F60 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.6)),
#       F80 = map_dbl(data, ~ get_F_at_spr_threshold(.x, 0.8)),
#       SSB    = map2_dbl(data, Fest,  ~ get_SSB_at_F(.x, .y)),
#       SSBMSY = map2_dbl(data, FMSY,  ~ get_SSB_at_F(.x, .y)),
#       SSB01  = map2_dbl(data, F01,   ~ get_SSB_at_F(.x, .y)),
#       SSB20  = map2_dbl(data, F20,   ~ get_SSB_at_F(.x, .y)),
#       SSB40  = map2_dbl(data, F40,   ~ get_SSB_at_F(.x, .y)),
#       SSB60  = map2_dbl(data, F60,   ~ get_SSB_at_F(.x, .y)),
#       SSB80  = map2_dbl(data, F80,   ~ get_SSB_at_F(.x, .y)),
#       FFMSY = Fest/FMSY, SSBSSBMSY = SSB/SSBMSY,
#       FF01 = Fest/F01, SSBSSB40 = SSB/SSB40) |>
#     dplyr::select(-data) |>
#     pivot_longer(-c(sim, spp), names_to = 'param', values_to = 'value')
# 
#   return(easi_stats)
# 
# }
# 
# # get_easi_stats_safe <- function(Bprod_summ, Fsusc_fshry) {
# #   # Compute F1 slope value
# #   
# #   
# #   # F1_value <- Bprod_summ |>
# #   #   dplyr::filter(x <= 0.02) |>
# #   #   group_by(sim, spp) |>
# #   #   mutate(slope = lag(((lead(ypr) - ypr) / (lead(x) - x)) * 0.1)) |>
# #   #   dplyr::filter(x == 0.02) |>
# #   #   dplyr::select(sim, spp, F1_slope = slope)
# #   
# #   # Get estimate of current fishing mortality
# #   totF <- Fsusc_fshry |>
# #     group_by(sim, spp) |>
# #     summarise(Fest = sum(Fage), susc = 1 - exp(-sum(Fage)), .groups = "drop")
# #   
# #   # Try one sim-spp at a time
# #   safe_results <- Bprod_summ |>
# #     group_by(sim, spp) |>
# #     tidyr::nest() |>
# #     #dplyr::left_join(F1_value, by = c("sim", "spp")) |>
# #     dplyr::left_join(totF, by = c("sim", "spp")) |>
# #     mutate(stats = purrr::pmap(
# #       list(data, Fest, FMSY = NA_real_, F01 = NA_real_),
# #       function(df, Fest, FMSY, F01) {
# #         tryCatch({
# #           FMSY <- df$x[which.max(df$ypr)]
# #           F01 <- df$x[which(df$slope <= F1_slope)[1]]
# #           
# #           Fsprs <- purrr::map_dbl(c(0.2, 0.4, 0.6, 0.8), function(thresh) {
# #             idx <- which.max(df$spr <= thresh)
# #             if (length(idx) == 0 || is.na(idx)) return(NA_real_)
# #             df$x[idx]
# #           })
# #           
# #           names(Fsprs) <- c("F20", "F40", "F60", "F80")
# #           
# #           get_SSB_at_F <- function(df, F_val) {
# #             if (is.na(F_val)) return(NA_real_)
# #             if (F_val == 0) return(0)
# #             idx <- which(df$x <= F_val)
# #             if (length(idx) == 0) return(NA_real_)
# #             df$SSBr[idx[which.max(df$x[idx])]]
# #           }
# #           
# #           SSB_vals <- purrr::map_dbl(
# #             c(Fest, FMSY, F01, Fsprs),
# #             ~ get_SSB_at_F(df, .x)
# #           )
# #           names(SSB_vals) <- c("SSB", "SSBMSY", "SSB01", "SSB20", "SSB40", "SSB60", "SSB80")
# #           
# #           FFMSY <- Fest / FMSY
# #           SSBSSBMSY <- SSB_vals["SSB"] / SSB_vals["SSBMSY"]
# #           
# #           tibble::tibble(
# #             param = c("FMSY", "F01", names(Fsprs), names(SSB_vals), "Fest", "FFMSY", "SSBSSBMSY"),
# #             value = c(FMSY, F01, Fsprs, SSB_vals, Fest, FFMSY, SSBSSBMSY)
# #           )
# #         }, error = function(e) {
# #           message("Skipping sim = ", unique(df$sim), ", spp = ", unique(df$spp), ": ", e$message)
# #           return(NULL)
# #         })
# #       }
# #     )) |>
# #     tidyr::unnest(cols = stats) |>
# #     dplyr::select(sim, spp, param, value)
# #   
# #   return(safe_results)
# # }
# 
# get_stat_idxs <- function(Bprod_summ, Fsusc_fshry) {
#   
#   Fest_vals <- Fsusc_fshry |>
#     group_by(sim, spp) |>
#     dplyr::summarise(Fest = sum(Fage, na.rm = T))
#   
#   idx_df <-
#     Bprod_summ |>
#     group_by(sim, spp) |>
#     nest() |>
#     left_join(Fest_vals, by = c("sim", "spp")) |>
#     mutate(FMSY = map_dbl(data, ~which.max(.x$ypr)),
#            F01 = map_dbl(data, ~which.min(abs(.x$slope - first(.x$slope) * 0.1))),
#            F20 = map_dbl(data, ~which.max(.x$spr <= 0.2)),
#            F40 = map_dbl(data, ~which.max(.x$spr <= 0.4)),
#            F60 = map_dbl(data, ~which.max(.x$spr <= 0.6)),
#            F80 = map_dbl(data, ~which.max(.x$spr <= 0.8)),
#            Fest = map2_dbl(data, Fest, ~which.min(abs(.x$x - .y)))) |>
#     dplyr::select(-data) |>
#     pivot_longer(-c(spp, sim), names_to = 'param', values_to = 'index')
#   
#   Fidx_df <- Bprod_summ |>
#     group_by(sim, spp) |>
#     mutate(rowid = row_number()) |>
#     ungroup() |>
#     right_join(idx_df, by = c('sim', 'spp')) |>
#     dplyr::filter(rowid == index) |>
#     dplyr::select(sim, spp, param, everything(), -rowid, -index)
#   
#   Fs_wide <-
#     Fidx_df |>
#     dplyr::select(sim, spp, param, x) |>
#     pivot_wider(names_from = param, values_from = x) 
#   
#   B_nested <-
#     Bprod_summ |>
#     group_by(sim, spp) |>
#     mutate(rowid = row_number()) |>
#     nest() |>
#     right_join(Fs_wide, by = c('sim', 'spp'))
#   
#   Bidx_df <-
#     B_nested |>
#     pivot_longer(cols = c(FMSY, F01, Fest, F20, F40, F60, F80),
#                  names_to = 'param', values_to = 'F_val') |>
#     mutate(param = recode(param, FMSY = "SSBMSY", F01 = "SSB01",
#                           Fest ="SSB", F20 = "SSB20", F40 = 'SSB40',
#                           F60 = "SSB60", F80 = "SSB80"),
#            index = map2_dbl(data, F_val, ~ which.min(abs(.x$x - .y)))) |>
#     dplyr::select(sim, spp, param, index, data) |>
#     unnest(data) |>
#     dplyr::filter(rowid == index) |>
#     dplyr::select(sim, spp, param, everything(), -rowid, -index) 
#   
#   idx_df <- bind_rows(Bidx_df, Fidx_df)
#   return(idx_df)
#   
# }
# 
# 
# #############################
# # `%||%` <- function(a, b) {
# #   if (is.null(a)) return(b)
# #   if (length(a) > 1) {
# #     b <- rep(b, length.out = length(a))
# #     a[is.na(a)] <- b[is.na(a)]
# #     return(a)
# #   }
# #   if (is.na(a)) return(b)
# #   return(a)
# # }
# 
# # vary_param <- function(x, name) {
# #   setting <- param_settings[[name]]
# #   if (is.null(setting)) return(x)  # No variation for this param
# # 
# #   vapply(x, function(val) {
# #     easi_param2(est = val,
# #                 dist = setting$dist %||% "fixed",
# #                 max = setting$max,
# #                 min = setting$min,
# #                 sd  = setting$sd,
# #                 se  = setting$se,
# #                 cv  = setting$cv)
# #   }, numeric(1))
# # }
# 
# prep_fishery_data <- function(fish_df, param_settings, verbose = TRUE) {
#   if (verbose) message("🎣 Preparing fishery simulations...")
#   
#   spp_list <- unique(fish_df$spp)
#   spp_data <- vector("list", length = length(spp_list))
#   names(spp_data) <- spp_list
#   
#   for (spp in spp_list) {
#     if (verbose) message("🟢 Processing species: ", spp)
#     
#     spp_df <- fish_df %>% filter(spp == !!spp)
#     fishery_data <- vector("list", length = nrow(spp_df))
#     
#     for (i in seq_len(nrow(spp_df))) {
#       row <- spp_df[i, ]
#       nsim <- row$nsim
#       
#       # Repeat row nsim times
#       sim_data <- row[rep(1, nsim), ] %>%
#         rowid_to_column("sim")
#       
#       # Apply vary_param to specified parameters
#       for (param in names(param_settings)) {
#         if (param %in% names(sim_data)) {
#           sim_data[[param]] <- vary_param2(sim_data[[param]], name = param, param_settings = param_settings)
#         }
#       }
#       
#       # Drop 'nsim' column
#       sim_data <- sim_data %>% ungroup |> select(-nsim)
#       
#       fishery_data[[i]] <- sim_data
#     }
#     
#     spp_data[[spp]] <- bind_rows(fishery_data)
#   }
#   
#   if (verbose) message("✅ Fishery simulations prepared.")
#   return(spp_data)
# }
# 
# prep_life_history <- function(lh_df, param_settings, verbose = TRUE) {
#   
#   
#   if (verbose) message("📥 Preparing life history simulations...")
#   
#   # Get unique length intervals
#   unique_len_ints <- lh_df %>%
#     pull(len_int) %>%
#     unique() %>%
#     na.omit() %>%
#     sort()
#   
#   if (verbose) message("🔢 Found length intervals: ", paste(unique_len_ints, collapse = ", "))
#   
#   # Precompute length bins for each len_int
#   bins_list <- map(unique_len_ints, function(li) {
#     max_Linf <- max(lh_df$Linf[lh_df$len_int == li], na.rm = TRUE)
#     seq(0, ceiling((max_Linf + li) / li) * li, by = li)
#   })
#   names(bins_list) <- as.character(unique_len_ints)
#   
#   if (verbose) message("📊 Length bins precomputed")
#   
#   spp_list <- lh_df$spp
#   spp_data <- vector("list", length = length(spp_list))
#   names(spp_data) <- spp_list
#   
#   for (i in seq_along(spp_list)) {
#     spp_row <- lh_df[i, ]
#     spp <- spp_row$spp
#     if (verbose) message("🟢 Processing species: ", spp)
#     
#     nsim <- as.integer(spp_row$nsim)
#     li <- spp_row$len_int
#     len_bins <- bins_list[[as.character(li)]]
#     
#     # Expand rows by simulation
#     sim_params <- spp_row[rep(1, nsim), ] %>%
#       rowid_to_column(var = 'sim')
#     
#     # Apply variation to parameters
#     for (param in names(param_settings)) {
#       if (param %in% names(sim_params)) {
#         sim_params[[param]] <- vary_param2(sim_params[[param]], name = param, param_settings = param_settings)
#       }
#     }
#     
#     # Ensure numeric and default values
#     sim_params <- sim_params %>%
#       mutate(
#         Linf = as.numeric(Linf),
#         L0 = as.numeric(L0 %||% 0),
#         len_int = li
#       )
#     
#     # Expand by length
#     expanded_df <- sim_params %>%
#       distinct(spp, sim, L0, Linf, .keep_all = TRUE) %>%
#       crossing(length = len_bins) %>%
#       #filter(length >= 0, length <= (Linf + li)) %>%
#       mutate(mn_length = length - (li / 2)) |>
#       #mutate(mn_length = (length + dplyr::lag(length, default = 0)) / 2) |>
#       mutate(spp = spp, max_age = max_age) |>
#       dplyr::filter(mn_length >= 0 & mn_length <= (Linf + (+ li / 2)))
#     
#     # Clean columns
#     expanded_df <- expanded_df %>%
#       select(-matches("_dist|_sd|_se|_min|_max")) %>%
#       select(spp, sim, length, mn_length, max_age, everything())
#     
#     spp_data[[spp]] <- expanded_df
#   }
#   
#   if (verbose) message("✅ All species prepared.")
#   return(spp_data)
# }
# 
# #param_settings <- list(
# #avail     = list(dist = "normal", sd = 0.05),
# #olap       = list(dist = "normal", sd = 0.1),
# #prm  = list(dist = "lognormal", cv = 0.2)
# #)
# 
# #x <- prep_life_history(lh_dat, param_settings = list())
# 
# prep_life_history2 <- function(lh_df, verbose = TRUE) {
#   if (verbose) message("📥 Preparing life history simulations...")
#   
#   # Identify variable parameters
#   varied_params <- lh_df %>%
#     dplyr::select(ends_with("_dist")) %>%
#     names() %>%
#     stringr::str_remove("_dist$")
#   
#   # Precompute length bins
#   unique_len_ints <- lh_df %>%
#     dplyr::pull(len_int) %>%
#     unique() %>%
#     na.omit() %>%
#     sort()
#   
#   if (verbose) message("🔢 Found length intervals: ", paste(unique_len_ints, collapse = ", "))
#   
#   bins_list <- purrr::map(unique_len_ints, function(li) {
#     max_Linf <- max(lh_df$Linf[lh_df$len_int == li], na.rm = TRUE)
#     seq(0, ceiling((max_Linf + li) / li) * li, by = li)
#   })
#   names(bins_list) <- as.character(unique_len_ints)
#   
#   if (verbose) message("📊 Length bins precomputed")
#   
#   spp_list <- lh_df$spp
#   spp_data <- vector("list", length = length(spp_list))
#   names(spp_data) <- spp_list
#   
#   for (i in seq_along(spp_list)) {
#     spp_row <- lh_df[i, ]
#     spp <- spp_row$spp
#     if (verbose) message("🟢 Processing species: ", spp)
#     
#     nsim <- as.integer(spp_row$nsim[[1]])
#     li <- spp_row$len_int
#     len_bins <- bins_list[[as.character(li)]]
#     
#     # Expand rows by simulation count
#     sim_params <- spp_row[rep(1, nsim), ] %>%
#       rowid_to_column(var = 'sim')
#     
#     # Apply variation via easi_param2
#     for (p in varied_params) {
#       base_val <- as.numeric(spp_row[[p]])
#       dist <- spp_row[[paste0(p, "_dist")]]
#       sd <- as.numeric(spp_row[[paste0(p, "_sd")]])
#       se <- as.numeric(spp_row[[paste0(p, "_se")]])
#       min <- as.numeric(spp_row[[paste0(p, "_min")]])
#       max <- as.numeric(spp_row[[paste0(p, "_max")]])
#       
#       sim_params[[p]] <- vapply(
#         seq_len(nsim),
#         function(j) {
#           if (j == 1) return(base_val)
#           easi_param2(est = base_val, dist = dist, max = max, min = min, sd = sd, se = se)
#         },
#         numeric(1)
#       )
#     }
#     
#     # Safely ensure numeric columns
#     sim_params <- sim_params %>%
#       dplyr::mutate(Linf = as.numeric(Linf),
#                     L0 = ifelse(is.na(L0), 0, as.numeric(L0)),
#                     len_int = li, max_age = max_age)
#     
#     # Expand by length bins
#     expanded_df <- sim_params %>%
#       dplyr::mutate(Linf = as.numeric(Linf), L0 = ifelse(is.na(L0), 0, as.numeric(L0)),
#                     len_int = li) |>
#       dplyr::distinct(spp, sim, L0, Linf, .keep_all = TRUE) %>%
#       tidyr::crossing(length = len_bins) %>%
#       dplyr::filter(length >= (L0 - li), length <= (Linf + li)) %>%
#       dplyr::mutate(mn_length = length + li / 2, spp = spp, max_age = max_age)
#     
#     # Clean up and reorder
#     expanded_df <- expanded_df %>%
#       dplyr::select(-dplyr::matches("_dist|_sd|_se|_min|_max")) %>%
#       dplyr::select(spp, sim, length, mn_length, max_age, dplyr::everything())
#     
#     spp_data[[spp]] <- expanded_df
#   }
#   
#   if (verbose) message("✅ All species prepared.")
#   return(spp_data)
# }
# 
# 
# get_easi_stats <- function(Bprod_summ, Fsusc_fshry) {
#   
#   # F1_value <-
#   #   Bprod_summ |>
#   #   dplyr::filter(x <= 0.02) |>
#   #   group_by(sim, spp) |>
#   #   mutate(slope = lag(((lead(ypr) - ypr)/(lead(x) - x))*0.1)) |>
#   #   dplyr::filter(x == 0.02) |>
#   #   dplyr::select(sim, spp, F1_slope = slope)
#   
#   F0_1_df <- Bprod_summ %>%
#     group_by(sim, spp) %>%
#     summarise(F1_slope = first(slope) * 0.1,
#       F01 = {idx <- which(slope <= F1_slope)
#         if (length(idx) == 0) NA_real_ else x[min(idx)]}, .groups = "drop")
#   
#   # Get est of F instantaneous (current F)
#   totF <-
#     Fsusc_fshry |>
#     group_by(sim, spp) |>
#     summarise(Fest = sum(Fage), susc = 1 - exp(-sum(Fage)))
#   #totF
#   
#   # Get FMSY estimate and cbind Finst to table
#   # ypr_stats <-
#   #   Bprod_summ |>
#   #   group_by(sim, spp) |>
#   #   nest() |>
#   #   left_join(F1_value, by = c('sim', 'spp')) |>
#   #   left_join(totF, by = c('sim', 'spp')) |>
#   #   mutate(FMSY = map(.x = data, .f = ~.$x[which.max(.$ypr)]), 
#   #          F10 = map(.x=data, .f = ~.$x[which.max(.$slope <= F1_slope)])) |>
#   #   dplyr::select(-c(data,susc)) |>
#   #   unnest()
#   
#   ypr_stats <-
#     Bprod_summ |>
#     left_join(F1_value, by = c("sim", "spp")) |>
#     group_by(sim, spp) |>
#     summarise(FMSY = x[which.max(ypr)],
#       F01 = x[which(slope <= F1_slope)[1]],  # first x where slope <= threshold
#       .groups = "drop") |>
#     left_join(totF, by = c('sim', 'spp'))
#   
#   stats <- 
#     Bprod_summ |>
#     group_by(sim, spp) |> 
#     nest() |>
#     left_join(ypr_stats, by = c('sim', 'spp')) |>
#     mutate(
#       # Calculate F values based on SPR thresholds
#       F20 = map_dbl(data, ~ .x$x[which.max(.x$spr <= 0.2)]),
#       F40 = map_dbl(data, ~ .x$x[which.max(.x$spr <= 0.4)]),
#       F60 = map_dbl(data, ~ .x$x[which.max(.x$spr <= 0.6)]),
#       F80 = map_dbl(data, ~ .x$x[which.max(.x$spr <= 0.8)]),
#       # Calculate SSB values at different F thresholds using helper function
#       SSB    = map2_dbl(data, Fest,  ~ get_SSB_at_F(.x, .y)),
#       SSBMSY = map2_dbl(data, FMSY,  ~ get_SSB_at_F(.x, .y)),
#       SSB01  = map2_dbl(data, F01,   ~ get_SSB_at_F(.x, .y)),
#       SSB20  = map2_dbl(data, F20,   ~ get_SSB_at_F(.x, .y)),
#       SSB40  = map2_dbl(data, F40,   ~ get_SSB_at_F(.x, .y)),
#       SSB60  = map2_dbl(data, F60,   ~ get_SSB_at_F(.x, .y)),
#       SSB80  = map2_dbl(data, F80,   ~ get_SSB_at_F(.x, .y)),
#       # Calculate ratio of F to FMSY
#       FFMSY = Fest / FMSY
#     ) |>
#     select(-data) |>  # Remove nested data list-column
#     unnest(cols = everything()) |>  # Flatten all columns
#     mutate(SSBSSBMSY = SSB / SSBMSY) |>  # Calculate ratio
#     pivot_longer(-c(sim, spp), names_to = 'param', values_to = 'value')  # Long format
#   
#   # stats <- 
#   #   Bprod_summ |>
#   #   group_by(sim, spp) |> 
#   #   nest() |>
#   #   left_join(ypr_stats, by = c('sim', 'spp')) |>
#   #   mutate(F20 = map(.x = data, .f = ~.$x[which.max(.$spr <= 0.2)]),
#   #          F40 = map(.x = data, .f = ~.$x[which.max(.$spr <= 0.4)]),
#   #          F60 = map(.x = data, .f = ~.$x[which.max(.$spr <= 0.6)]),
#   #          F80 = map(.x = data, .f = ~.$x[which.max(.$spr <= 0.8)]),
#   #          SSB = map(.x = data, .f = ~.$SSBr[which.min(.$x <= Fest)]),
#   #          SSBMSY = map(.x = data, .f = ~.$SSBr[which.min(.$x <= FMSY)]),
#   #          SSB01 = map(.x = data, .f = ~.$SSBr[which.min(.$x <= F01)]),
#   #          SSB20 = map(.x = data, .f = ~.$SSBr[which.min(.$x <= F20)]),
#   #          SSSB40 = map(.x = data, .f = ~.$SSBr[which.min(.$x <= F40)]),
#   #          SSB60 = map(.x = data, .f = ~.$SSBr[which.min(.$x <= F60)]),
#   #          SSB80 = map(.x = data, .f = ~.$SSBr[which.min(.$x <= F80)]),
#   #          FFMSY = Fest/FMSY) |>
#   #   dplyr::select(-data) |>
#   #   unnest() |>
#   #   mutate(SSBSSBMSY = SSB/SSBMSY) |>
#   #   pivot_longer(-c(sim,spp), names_to = 'param', values_to = 'value')
#   return(stats)
# }
# 
# safe_F_threshold <- function(df, threshold, fallback = NA) {
#   filtered <- df[df$spr <= threshold, ]
#   if (nrow(filtered) == 0) return(NA)
#   min(filtered$x)
# }
# 
# get_SSB_at_F <- function(df, F_val) {
#   idx <- which(df$x <= F_val)
#   if(length(idx) == 0) return(NA_real_)
#   df$SSBr[idx[which.max(df$x[idx])]]
# }
# 
# get_ypr_spr <- function(Fsusc_len, lh_ins, Fmax = 10, dF = 0.02) {
#   df <- left_join(Fsusc_len, lh_ins, by = c('spp', 'sim', 'mn_length'))
#   
#   # Fishing mortality values
#   F_vals <- seq(0, Fmax, by = dF)
#   nF <- length(F_vals)
#   
#   # Sort input by mn_length
#   df <- df[order(df$mn_length), ]
#   
#   # Extract vectors
#   Fsusc <- df$Fsusc
#   Fsusc_min <- unique(df$Fsusc_min)
#   M <- df$M
#   dT <- df$dT
#   wt <- df$wt
#   mat <- df$mat
#   len <- df$length
#   mn_length  <- df$mn_length
#   
#   # Output storage
#   ypr_out   <- numeric(nF)
#   ssbr_out  <- numeric(nF)
#   
#   for (i in seq_along(F_vals)) {
#     F_now <- F_vals[i]
#     
#     z <- (Fsusc * F_now + M) * dT
#     
#     # YPR calculation (length-based)
#     ypr_term <- z
#     yprF <- ifelse(
#       ypr_term == 0, 0,
#       (wt * (Fsusc * F_now) / (Fsusc * F_now + M)) *
#         ((1 - exp(-ypr_term)) * exp(-cumsum(ypr_term)))
#     )
#     ypr_out[i] <- mean(yprF, na.rm = TRUE)
#     
#     # SSB per recruit
#     SSBlc <- rep(NA_real_, length(z))
#     SSBlc[Fsusc == 0] <- 0
#     SSBlc[len == Fsusc_min] <- 1
#     cond <- which(Fsusc > 0 & mn_length > Fsusc_min)
#     if (length(cond) > 0) {
#       SSBlc[cond] <- exp(-dplyr::lag(ypr_term, default = NA_real_)[cond])
#     }
#     SSBlc_cum <- cumprod(replace_na(SSBlc, 1))
#     SSBr <- ifelse(Fsusc == 0, NA_real_, wt * mat * SSBlc_cum)
#     ssbr_out[i] <- sum(SSBr, na.rm = TRUE)
#   }
#   
#   # Combine results
#   sim_id <- df$sim[1]
#   spp_id <- df$spp[1]
#   ssbr0  <- ssbr_out[which(F_vals == 0)][1]  # SSBr at F = 0
#   
#   result <- tibble(
#     sim    = sim_id,
#     spp    = spp_id,
#     x      = F_vals,
#     ypr    = ypr_out,
#     SSBr   = ssbr_out,
#     slope  = (lead(ypr_out) - ypr_out) / (lead(F_vals) - F_vals),
#     spr    = ssbr_out / ssbr0
#   )
#   
#   return(result)
# }
# 
# # ypr_results <- tmp2 %>%
# #   group_by(sim, spp) %>%
# #   group_split() %>%
# #   map_dfr(get_ypr_spr, Fmax = 10, dF = 0.02)
# 
# 
# # Helper to apply easi_param2 vectorized over a column
# # vary_param <- function(x, name, param_settings) {
# #   setting <- param_settings[[name]]
# #   
# #   if (is.null(setting)) {
# #     return(x)  # No variation
# #   }
# #   
# #   x <- as.numeric(x) 
# #   
# #   vapply(x, function(val) {
# #     easi_param2(est = val,
# #                 dist = setting$dist %||% "fixed",
# #                 max = setting$max,
# #                 min = setting$min,
# #                 sd  = setting$sd,
# #                 se  = setting$se,
# #                 cv  = setting$cv)
# #   }, numeric(1))
# # }
# 
# rnorm_clipped <- function(mean, sd, min = 0, max = 1) {
#   val <- rnorm(1, mean, sd)
#   val <- max(min(val, max), min)  # Clip to bounds
#   return(val)
# }
# 
# vary_param2 <- function(x, name, param_settings) {
#   setting <- param_settings[[name]]
#   
#   if (is.null(setting)) return(x)
#   
#   x <- as.numeric(x)
#   
#   vapply(x, function(val) {
#     # Evaluate min/max if they are functions of val
#     min_val <- if (is.function(setting$min)) setting$min(val) else setting$min
#     max_val <- if (is.function(setting$max)) setting$max(val) else setting$max
#     
#     easi_param2(
#       est  = val,
#       dist = setting$dist %||% "fixed",
#       min  = min_val,
#       max  = max_val,
#       sd   = setting$sd,
#       se   = setting$se,
#       cv   = setting$cv
#     )
#   }, numeric(1))
# }
# 
# # Helper to apply easi_param2 vectorized over a column
# # vary_param <- function(x, name, param_settings) {
# #   setting <- param_settings[[name]]
# #   
# #   if (is.null(setting)) {
# #     return(x)  # No variation
# #   }
# #   
# #   x <- as.numeric(x) 
# #   
# #   vapply(x, function(val) {
# #     easi_param2(est = val,
# #                 dist = setting$dist %||% "fixed",
# #                 max = setting$max,
# #                 min = setting$min,
# #                 sd  = setting$sd,
# #                 se  = setting$se,
# #                 cv  = setting$cv)
# #   }, numeric(1))
# # }
# 
# 
# easi_param2 <- function(est, dist = 'fixed', max = NA, min = NA, sd = NA, se = NA, cv = NA) {
#   library(truncnorm)
#   if (is.na(est) || est == "") return(NA)
#   
#   if (dist == 'fixed') {
#     z <- est
#     
#   } else if (dist == 'normal') {
#     z <- rnorm(1, mean = est, sd = sd)
#     
#   } else if (dist == 'normal01') {
#     z <- rnorm_clipped(est, sd, min = 0, max = 1)
#     
#   } else if (dist == 'uniform') {
#     z <- runif(1, min = min, max = max)
#     
#   } else if (dist == 'lognormal') {
#     if (is.na(cv)) stop("Must provide cv for lognormal")
#     sdlog <- sqrt(log(1 + cv^2))
#     meanlog <- log(est) - 0.5 * sdlog^2
#     z <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
#     
#   } else if (dist == 'truncnorm') {
#     if(is.na(sd)) stop("missing SD for truncnorm")
#     z <- rtruncnorm(1, a=min, b=max, mean=est, sd=sd)
#   } 
#   else {
#     z <- NA
#   }
#   
#   return(z)
# }
# 
# # easi_tidy function ---- 
# # This function takes a df with first column param, and then following column to be spp or fleets
# # and puts them into a nested list of length spp/fleet and length(params) within each element
# # Used in fleet_data script
# easi_tidy <- function(df) {
#   params <- df$param
#   nms <- names(df)[-1]  # Exclude the 'param' column
#   y <- list()
#   
#   map(nms, ~ {y[[.x]] <<- setNames(as.list(df[[.x]]), params)})
#   
#   return(y)
# }
# 
# # easi_param function ----
# # This function takes in lh and fsry params and applies a distribution to them (fixed, normal, or uniform)
# # Each param has a row for: estimate, distribution, min, max, sd, and se from input csv.
# easi_param <- function(x, fsry = NA, pmtr = NA) { 
#   if(is.na(fsry) || length(fsry) == 0) {
#     y <- x
#   } else {
#     # If fsry is a character vector with more than one element,
#     # select only those columns of x that match the elements of fsry
#     if(length(fsry) > 1) {
#       y <- x[fsry]
#     } else {
#       # If fsry is a single character string, select the corresponding column of x
#       y <- x[[fsry]]
#     }
#   }
#   nms <- c(pmtr, paste0(pmtr,'_dist'), paste0(pmtr,'_max'), paste0(pmtr,'_min'), paste0(pmtr,'_sd'), paste0(pmtr,'_se'))
#   y <- y[nms] |> map(type.convert, as.is = TRUE)
# 
#   y <- y[order(names(y))]
#   p = y[[1]]; dist = y[[2]]; max = y[[3]]; min = y[[4]]; sd = y[[5]]; if(length(y) == 6){se = y[[6]]}
#   if(is.na(p) | p == ""){return(NA)}
#   if(dist == 'fixed') {z <- p}
#   if(dist == 'normal') {z <- rnorm(1,p,sd)}
#   if(dist == 'uniform') {z <- runif(1,min,max)}
#   return(as.numeric(z))
# }
# 
# # easi_param2 fuinction ----
# # same as above recoded to diy within mutate etc
# # This function takes in lh and fsry params and applies a distribution to them (fixed, normal, or uniform)
# # Each param has a row for: estimate, distribution, min, max, sd, and se from input csv.
# 
# # easi_param2 <- function(est, dist = 'fixed', max = NA, min = NA, sd = NA, se = NA) { 
# #   # if(is.na(fsry) || length(fsry) == 0) {
# #   #   y <- x
# #   # } else {
# #   #   # If fsry is a character vector with more than one element,
# #   #   # select only those columns of x that match the elements of fsry
# #   #   if(length(fsry) > 1) {
# #   #     y <- x[fsry]
# #   #   } else {
# #   #     # If fsry is a single character string, select the corresponding column of x
# #   #     y <- x[[fsry]]
# #   #   }
# #   # }
# #   # nms <- c(pmtr, paste0(pmtr,'_dist'), paste0(pmtr,'_max'), paste0(pmtr,'_min'), paste0(pmtr,'_sd'), paste0(pmtr,'_se'))
# #   # y <- y[nms] |> map(type.convert, as.is = TRUE)
# #   # 
# #   # y <- y[order(names(y))]
# #   # p = y[[1]]; dist = y[[2]]; max = y[[3]]; min = y[[4]]; sd = y[[5]]; if(length(y) == 6){se = y[[6]]}
# #   if(is.na(est) | est == ""){return(NA)}
# #   if(dist == 'fixed') {z <- est}
# #   if(dist == 'normal') {z <- rnorm(1,est,sd)}
# #   if(dist == 'uniform') {z <- runif(1,min,max)}
# #   return(as.numeric(z))
# # }
# # 
# # easi_param2 <- function(est, dist = 'fixed', max = NA, min = NA, sd = NA, se = NA, cv = NA) {
# #   if (is.na(est) || est == "") return(NA)
# #   
# #   if (dist == 'fixed') {
# #     z <- est
# #   } else if (dist == 'normal') {
# #     z <- rnorm(1, mean = est, sd = sd)
# #   } else if (dist == 'uniform') {
# #     z <- runif(1, min = min, max = max)
# #   } else if (dist == 'lognormal') {
# #     if (is.na(cv)) stop("Must provide cv for lognormal")
# #     sdlog <- sqrt(log(1 + cv^2))
# #     meanlog <- log(est) - 0.5 * sdlog^2
# #     z <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
# #   } else {
# #     z <- NA
# #   }
# #   
# #   return(z)
# # }
# 
# # Life history functions ----
# 
# # easi_growth function ----
# # easi_growth converts age from length across various equations
# # easi_growth - calc age based on various growth eqs and associated params
# easi_growth2 <- function(model, len, Linf, K, t0 = NA, L0 = NA, D = NA, m = NA,
#                         age1 = NA, age2 = NA, len1 = NA, len2 = NA,
#                         a = NA, b = NA, alpha = NA, g = NA) {
#   
#   n <- length(len)
#   
#   # Helper to recycle scalars to vector length
#   r <- function(x) if (length(x) == 1) rep(x, n) else x
#   
#   model <- r(model); Linf <- r(Linf); K <- r(K); t0 <- r(t0); L0 <- r(L0)
#   D <- r(D); m <- r(m); age1 <- r(age1); age2 <- r(age2)
#   len1 <- r(len1); len2 <- r(len2); a <- r(a); b <- r(b); alpha <- r(alpha); g <- r(g)
#   
#   out <- numeric(n)
#   
#   for (i in seq_len(n)) {
#     out[i] <- switch(model[i],
#                      "SVBGF"     = log(1 - (len[i]/Linf[i])) / -K[i] + t0[i],
#                      "GVBGF"     = log(1 - (len[i]/Linf[i])^D[i]) / -K[i] + t0[i],
#                      "VBL0"      = 1 / K[i] * log((Linf[i] - L0[i]) / (Linf[i] - len[i])),
#                      "Gompertz"  = -log(-log(len[i] / Linf[i]) / (1/K[i])) / K[i] + t0[i],
#                      "Gompertz2" = -1 / K[i] * log(1 - log(len[i] / L0[i]) / log(Linf[i] / L0[i])),
#                      "Richards"  = -log(((Linf[i]/len[i])^m[i] - 1) / t0[i]) / K[i],
#                      "Schnute"   = age1[i] - (1/K[i]) * log((len[i] - Linf[i]) / (len1[i] - Linf[i])),
#                      "Schnute2"  = age1[i] + log(1 - ((len[i]^b[i] - len1[i]^b[i]) /
#                                                         (len2[i]^b[i] - len1[i]^b[i])) *
#                                                    (1 - exp(-a[i] * (age2[i] - age1[i])))) / (-a[i]),
#                      "Logistic"  = (1 / K[i]) * log((Linf[i] - alpha[i]) / (Linf[i] - len[i])),
#                      "Logistic2" = (1 / g[i]) * log((len[i] * L0[i] - len[i] * Linf[i]) /
#                                                       (len[i] * L0[i] - Linf[i] * L0[i])),
#                      "Logistic3" = (1 / K[i]) * log((len[i] * L0[i] - len[i] * Linf[i]) /
#                                                       (len[i] * L0[i] - Linf[i] * L0[i])),
#                      "Robertson" = (t0[i] - log((Linf[i] / len[i]) - 1)) / K[i],
#                      NA_real_
#     )
#   }
#   
#   return(out)
# }
# 
# # age1 and len1 for ALV only - update!
# 
# 
# easi_growth <- function(model = 'SVBGF', len = NULL, Linf = NULL, K = NULL,
#                         t0 = NULL, L0 = NULL, D = NULL, m = NULL, age1 = NULL, age2 = NULL, 
#                         len1 = NULL, len2 = NULL, a = NULL, b = NULL, alpha = NULL, g = NULL) {
#   #if(is.na(model)) {model <- 'SVGBF'}
#   if(model == 'SVBGF') { age <- log(1 - (len/Linf)) / -K + t0} # len, Linf, K, t0
#   if(model == 'GVBGF') { age <- log(1 - (len/Linf)^D)  / -K + t0} # len, Linf, D, K, t0
#   if(model == 'VBL0') { age <- 1 / K * log((Linf - L0) / (Linf - len))} #K, Linf, L0, len
#   #if(model == 'VBL0') { age <- 1 / K * log((Linf - len) / (Linf - L0))} #K, Linf, L0, len
#   if(model == 'Gompertz') { age <- -log(-log(len / Linf) / (1/K)) / K + t0} # len, Linf, K, t0
#   if(model == 'Gompertz2') {age <- -1 / K * log(1 - log(len / L0) / log(Linf / L0))}
#   if(model == 'Richards') { age <- -log(((Linf/len)^m-1) / t0) / K} # len, Linf, m, t0, K
#   if(model == 'Schnute') {age <- age1 - (1/K) * log((len - Linf) / (len1 - Linf))}
#   if(model == 'Schnute2') {age <- age1 + log(1 - ((len^b - len1^b) / (len2^b - len1^b)) * (1 - exp(-a * (age2 - age1)))) / (-a)}
#   if(model == 'Logistic') { age <- (1 / K) * log((Linf - alpha) / (Linf - len))}
#   if(model == 'Logistic2') {age <- (1 / g) * log((len * L0 - len * Linf) / (len * L0 - Linf * L0))} # silky                           
#   if(model == 'Logistic3') { age <- (1 / K) * log((len * L0 - len * Linf) / (len * L0 - Linf * L0))}#c plumbeus curve
#   if(model == 'Robertson') {age <- (t0 - log((Linf / len) - 1)) / K} # dusky = b param input as t0 here.
#   #logistic eq from Geraghty et al 2014 - len <- (Linf * L0 * exp(K * age))/(Linf + L0*(exp(K*age)-1))
#   #logistic from grant et al 2018 silky - len <- (Linf * L0 *exp(g*age)) / (Linf + L0 * (exp(g*age)-1))                                 
#   return(age)
# }
# 
# # tmp <- data.frame(len = seq(from = 0, to = 223, by = 4)) |>
# #   mutate(age = easi_growth(model = 'SVBGF', t0 = -4.5, Linf = 223, K = 0.1, len = len),
# #          dT = deltaT(K = 0.1, Linf = 223, len1 = len, len_int = 4),
# #          dT2 = deltaT2(K=0.1, Linf = 223, len = len, len_int = 4),
# #          dT3 = deltaT3(K = 0.1, Linf = 223, len1 = len, len_int = 4))
# # 
# # tmp |>
# #   ggplot() +
# #   geom_point(aes(age, len))
# 
# # easi_growth(model = 'VBL0', len = 1000, Linf = 1582, K = 0.021, t0 = 64,
# #             L0 = 64, D = NA, m = NA)
# 
# # length-weight function ----
# # Calculate weight from length. Beware units of length/weight from paper.
# lw <- function(a,len,b) {
#   wt <- (a*len^b)#/1000
#   return(wt)
# }
# 
# # data.frame(len = 0:308 ) |>
# #   mutate(wt = lw(a = 0.0000104, len = len, b = 2.9)) |>
# #   ggplot() +
# #   aes(len, wt) +
# #   geom_line()
# 
# # Delta-T function ----
# # Calculate deltaT - time it takes to grow from one len/age bin to the next
# # As written in easifish excel sheet, slightly varies from deltaT equation
# # deltaT <- function(K,Linf, len1, len_int) {
# #   len2 <- len1 - (len_int)
# #   
# #   dT <- ifelse(len1 == 0 | len1 > Linf, 0, 
# #                (1/K) * log((Linf - len2)/(Linf - len2 - (len1 - len2))))
# #   return(dT)
# # }
# 
# deltaT <- function(K, Linf, len1, len_int) {
#   len2 <- len1 - len_int
#   
#   # Cap len1 just below Linf to avoid instability
#   len1 <- ifelse(len1 >= Linf, Linf - 0.01, len1)
#   
#   dT <- ifelse(len1 <= 0 | len2 <= 0, 0, 
#                (1/K) * log((Linf - len2)/(Linf - len1)))
#   return(dT)
# }
# 
# # Delta-T2 fucntion ----
# # Written as in Griffiths et al 2023 report equation and chen and gordon
# # Above assigns delta T of prev len bin to current, this one doesnt. ie one row offset.
# # deltaT2 <- function(K,Linf, len, len_int) {
# #   len2 <- len - len_int
# #   dT <- ifelse(len == 0 | len > Linf, Linf-0.1,
# #                (1/K) * log((Linf - len2)/(Linf - len2 - len_int)))
# #   return(dT)
# # }
# 
# deltaT2 <- function(K, Linf, len, len_int) {
#   # Apply rules
#   len_adj <- case_when(
#     len == 0      ~ NA_real_,          # handle 0 separately below
#     len > Linf    ~ Linf - 0.1,
#     TRUE          ~ len)
#   
#   len2 <- len_adj - len_int
#   
#   dT <- ifelse(len == 0,0,
#                (1 / K) * log((Linf - len2) / (Linf - len2 - len_int)))
#   
#   return(dT)
# }
# 
# 
# 
# # Maturity function ----
# # Calculate maturity using L50, r and/or Lm 
# # Setup so that r can't be negative
# 
# easi_mat <- function(len, L50 = NA, r = NA, Lm = NA) {
#   if(!is.na(L50) & !is.na(r)) {mat <- 1 / (1 + exp(-r * (len - L50)))}
#   if(!is.na(L50) & is.na(r)) { mat <- ifelse(L50 >= len, 0, 1)}
#   if(is.na(L50) & is.na(r) & !is.na(Lm)) {mat <- ifelse(Lm > len, 0, 1)}
#   return(mat)
# }
# 
# #
# easi_mat2 <- function(len, L50 = NA, r = NA, Lm = NA) {
#   if (!all(is.na(L50)) && !all(is.na(r))) {
#     # Logistic curve
#     return(1 / (1 + exp(-r * (len - L50))))
#   }
#   
#   if (!all(is.na(L50)) && all(is.na(r))) {
#     # Knife-edge at L50
#     return(ifelse(len < L50, 0, 1))
#   }
#   
#   if (all(is.na(L50)) && all(is.na(r)) && !all(is.na(Lm))) {
#     # Knife-edge at Lm
#     return(ifelse(len < Lm, 0, 1))
#   }
#   
#   # Default: unknown params
#   return(rep(NA_real_, length(len)))
# }
# 
# # x <- data.frame(length = 0:320) |>
# #   mutate(mat1 = easi_mat(len = length, L50 = 285.3, r = 0.23))
# 
# # length 2 length function ----
# # Function that does length-length conversions for various shark species
# # Converts from total length (TL), fork length (UF), disk width (DW) and precaudal length (PC)
# # where available.
# # This function converts from an input length type to desired output length type.
# # If 2 of 3 functions are available, this function will do a double conversion to get to the required length type.
# 
# length_2_length <- function(spp, value, input_type, output_type) {
#   # Vectorize the whole function for convenience
#   # Convert vectors element-wise
#   mapply(function(sp, val, in_type, out_type) {
#     if (is.na(in_type) || is.na(out_type)) return(NA)
#     if (in_type == "DW") return(val)
#     if (in_type == out_type) return(val)
#     
#     direct_conv <- function(from, to, spp, val) {
#       if (from == "TL" && to == "UF") return(total_2_fork(spp, val))
#       if (from == "TL" && to == "PC") return(total_2_precaudal(spp, val))
#       if (from == "UF" && to == "TL") return(fork_2_total(spp, val))
#       if (from == "UF" && to == "PC") return(fork_2_precaudal(spp, val))
#       if (from == "PC" && to == "TL") return(precaudal_2_total(spp, val))
#       if (from == "PC" && to == "UF") return(precaudal_2_fork(spp, val))
#       return(NA)
#     }
#     
#     # Try direct conversion first
#     res <- direct_conv(in_type, out_type, sp, val)
#     if (!is.na(res)) return(res)
#     
#     # Try intermediates conversion
#     length_types <- c("TL", "UF", "PC")
#     intermediates <- setdiff(length_types, c(in_type, out_type))
#     
#     for (int_type in intermediates) {
#       step1 <- direct_conv(in_type, int_type, sp, val)
#       if (!is.na(step1)) {
#         step2 <- direct_conv(int_type, out_type, sp, step1)
#         if (!is.na(step2)) return(step2)
#       }
#     }
#     
#     # If still no result, return NA
#     return(NA)
#   }, spp, value, input_type, output_type)
# }
# 
# # total 2 fork function ----
# # Function to convert total to fork length for various species
# # Need FAO species code and total length measure
# # CCL, FAL, DUS, BLR
# 
# # data.frame(UF = 0:250) |>
# #   mutate(DUS = fork_2_total('DUS', UF),
# #          FAL = fork_2_total('FAL', UF),
# #          BLR = fork_2_total('BLR', UF),
# #          CCL = fork_2_total('CCL', UF)) |>
# #   pivot_longer(-UF, names_to = 'sp', values_to = 'TL') |>
# #   ggplot() +
# #   aes(UF, TL, col = sp) +
# #   geom_line()
# 
# total_2_fork <- function(spp, TL) {
#   UF <- NA
#   #ALS
#   if (spp == 'ALV') { UF <- 0.5474 * TL + 7.0262}  # Alopias vulpinus, Kohler et al. 1995
#   if (spp == 'AML') { UF <- (TL - 7.67)/1.13 } # Carcharhinus amblyrhyncos, Bradley et al. 2017
#   if (spp == 'BLR') { UF <- (TL - 49.06) / 1.15}  # Carcharhinus melanopterus, Chin et al. 2013
#   if (spp == 'BRO') { UF <- (TL - 22.544) / 0.83}  # Carcharhinus brachyurus, Drew et al. 2016
#   if (spp == 'BSH') { UF <- (0.838 * TL) - 1.615}  # Prionace glauca,  Francis & Duffy 2005
#   if (spp == 'BTH') { UF <- (TL + 7.2529) / 1.7273}  # Alopias superciliosus, Calle-Moran et al. 2023
#   if (spp == 'CCA') { UF <- 0.8074 * TL + 7.7694}  # Carcharhinus altimus, Kohler et al. 1995
#   #CCE
#   #CCG
#   if (spp == 'CCL') { UF <- (TL - 0.913) / 1.235}  # Carcharhinus limbatus, Killam & Parsons 1989
#   if (spp == 'CCP') { UF <- (TL - 2.747) / 1.206}  # Carcharhinus plumbeus, Geraghty et al. 2014
#   #CYW
#   if (spp == 'DUS') { UF <- (TL - 4.226) / 1.203}  # Carcharhinus obscurus, Geraghty et al. 2014
#   if (spp %in% c('FAL', 'CCG', 'ALS')) { UF <- (TL - 2.36) / 1.21}   # Carcharhinus falciformis, Joung et al 2008
#   #ISB
#   #LMA
#   #LMD
#   if (spp == 'OCS') { UF <- 0.817 * TL - 1.875}  # Carcharhinus longimanus, Joung et al. 2016
#   #PLS*
#   if (spp == 'POR') { UF <- 0.8896 * TL + 0.3369}  # Lamna nasus, Francis 2013
#   if (spp == 'PSK') { UF <- (TL - 5.168) / 1.1104}  # Pseudocarcharias kamoharai, Calle-Moran et al. 2024
#   if (spp == 'PTH') { UF <- (TL - 123.12) / 1.85}   # Alopias pelagicus, White 2007
#   if (spp == 'RHN') { UF <- (TL - 26.481) / 1.063 }   # Rhincodon typus, Wintner 2000
#   #RMB*
#   if (spp == 'RHU') { UF <- (TL - 4.46679) / 1.1538}  # Rhizoprionodon longurio, Marquez-Farias et al. 2005
#   #SCK
#   if (spp == 'SMA') { UF <- (TL * 0.9268) - 1.7101}   # Isurus oxyryinchus, Kohler et al. 1995
#   if (spp == 'SPK') { UF <- (TL - 3.58) / 1.29}  # Sphyrna mokarran, Stevens & Lyle 1989
#   if (spp == 'SPL') { UF <- (TL - 1.28) / 1.3}  # Sphyrna lewini, Stevens & Lyle 1989
#   if (spp == 'SPZ') { UF <- 0.817 * TL - 7.0834}  # Sphyrna zygaena, Chow 2004
#   if (spp == 'SSQ') { UF <- (TL - 2.74) / 1.08}  # Zameus squamulosus, Irvine et al. 2006
#   if (spp == 'TIG') { UF <- (TL - 22.607) / 1.096}  # Galeocerdo cuvier, Holmes et al. 2015
#   if (spp == 'WSH') { UF <- 0.94 * TL - 5.74}  # Carcharodon carcharias, Natanson & Skomal 2015
#   return(UF)
# }
# 
# # fork 2 total function ----
# # Function to convert fork to total length for various species
# # Need FAO species code and fork length measure
# fork_2_total <- function(spp, UF) {
#   TL <- NA
#   if (spp == 'ALV') { TL <- (UF - 7.0262) / 0.5474 } # Alopias vulpinus, Kohler et al. 1995
#   if (spp == 'AML') { TL <- 7.67 + 1.13 * UF }  # Carcharhinus amblyrhyncos, Bradley et al. 2017
#   if (spp == 'BLR') { TL <- 1.15 * UF + 49.06 } # Carcharhinus melanopterus, Chin et al. 2013
#   if (spp == 'BRO') { TL <- 0.83 * UF + 22.544 }# Carcharhinus brachyurus, Drew et al. 2016
#   if (spp == 'BSH') { TL <- (1.615 + UF) / 0.838}  # Prionace glauca,  Francis & Duffy 2005
#   if (spp == 'BTH') { TL <- 1.7273 * UF - 7.2529 }    # Alopias superciliosus, Calle-Moran et al. 2023
#   if (spp == 'CCA') { TL <- (UF - 7.7694) / 0.8074 }  # Carcharhinus altimus, Kohler et al. 1995
#   if (spp == 'CCL') { TL <- 1.235 * UF + 0.913 }# Carcharhinus limbatus, Killam & Parsons 1989
#   if (spp == 'CCP') { TL <- 1.206 * UF + 2.747 }# Carcharhinus plumbeus, Geraghty et al. 2014
#   if (spp == 'CYW') { TL <- 1.1 * UF } # guess
#   if (spp == 'DUS') { TL <- 1.203 * UF + 4.226 }# Carcharhinus obscurus, Geraghty et al. 2014
#   if (spp %in% c('FAL', 'CCG', 'ALS')) { TL <-  1.21 * UF + 2.36 } # Carcharhinus falciformis, Joung et al. 2008
#   if (spp == 'OCS') { TL <- (UF + 1.875) / 0.817 }  # Carcharhinus longimanus, Joung et al. 2016
#   if (spp == 'POR') { TL <- (UF - 0.3369) / 0.8896 } # Lamna nasus, Francis 2013
#   if (spp == 'PSK') { TL <- 1.1104 * UF + 5.168 }    # Pseudocarcharias kamoharai, Calle-Moran et al. 2024
#   if (spp == 'PTH') { TL <- (1.85 * UF) + 123.12 }  # Alopias pelagicus, White 2007
#   if (spp == 'RHN') { TL <- (UF * 1.063) + 26.481}   # Rhincodon typus, Wintner 2000
#   if (spp == 'RHU') { TL <- 1.1538 * UF + 4.46679 }  # Rhizoprionodon longurio, Marquez-Farias et al. 2005
#   if (spp == 'SMA') { TL <- (UF + 1.7101) / 0.9286} # Isurus oxyryinchus, Kohler et al. 1995
#   if (spp == 'SPK') { TL <- 1.29 * UF + 3.58 } # Sphyrna mokarran, Stevens & Lyle 1989
#   if (spp == 'SPL') { TL <- 1.3 * UF + 1.28 }  # Sphyrna lewini, Stevens & Lyle 1989
#   if (spp == 'SPZ') { TL <- (UF + 7.0834) / 0.817 }  # Sphyrna zygaena, Chow 2004
#   if (spp == 'SSQ') { TL <- 1.08 * UF + 2.74 } # Zameus squamulosus, Irvine et al. 2006
#   if (spp == 'TIG') { TL <- 1.096 * UF + 22.607 } # Galeocerdo cuvier, Holmes et al. 2015
#   if (spp == 'WSH') { TL <- (UF + 5.74) / 0.94 }# Carcharodon carcharias, Natanson & Skomal 2015
#   
#   return(TL)
# }
# 
# # precaudal 2 total function ----
# # Function to convert precaudal to total length for various species
# # Need FAO species code and precaudal length measure
# precaudal_2_total <- function(spp, PC) {
#   TL <- NA
#   if(spp == 'AML') { TL <- 17.88 + 1.17 * PC} # Carcharhinus amblyrhynchos, Bradley et al. 2017
#   if(spp == 'BRO') { TL <- (1.364 * PC) - 35.924} # Carcharhinus brachyurus, Drew et al. 2016
#   if(spp == 'BSH') { TL <- (3.75 + PC) / 0.78} # Prionace glauca,  Fujinami et al. 2017
#   if(spp == 'BTH') { TL <- (1.8232 * PC) + 0.7085} # Alopias superciliosus, Calle-Moran et al. 2023
#   if(spp == 'CCE') { TL <- (0.916 + PC) / 0.81} # Carcharhinus leucas, Cliff & Dufley, 1991
#   if(spp == 'CCG') { TL <- (1.17034 + PC) / 0.75865} # Carcharhinus galapagensis, Wetherbee et al. 1996
#   if(spp == 'CCP') { TL <- (1.316 * PC) + 4.566} # Carcharhinus plumbeus, Geraghty et al. 2014
#   if(spp == 'DUS') { TL <- (1.305 * PC) + 8.021} # Carcharhinus obscurus, Geraghty et al. 2014
#   if(spp == 'FAL') { TL <-  (1.31* PC) + 3.64} # Carcharhinus falciformis, Joung et al. 2008
#   if(spp == 'LMA') { TL <- (2.13 + PC) / 0.84} # Isurus paucus, Semba et al. 2009
#   if(spp == 'LMD') { TL <- (1.15 * PC) + 15.19} # Lamna ditropis, Goldman & Musick, 2006
#   if(spp == 'OCS') { TL <- (6.019 + PC) / 0.755} # Carcharhinus longimanus, Joung et al. 2016
#   if(spp == 'PTH') { TL <- 2.34 + (1.93 * PC)} # Alopias pelagicus, Liu et al. 1999
#   if(spp == 'RHN') { TL <- (1.148 * PC) + 0.262} # Rhincodon typus, Hsu et al. 2012
#   if(spp == 'SMA') { TL <- (2.13 + PC) / 0.84} # Isurus oxyrinchus, Semba et al. 2009
#   if(spp == 'SPZ') { TL <- 4.3873 + 1.3531 * PC} # Sphyrna zygaena, Lopez Martinez et al. 2020
#   if(spp == 'TIG') { TL <- 34.321 + 1.159 * PC} # Galeocerdo cuvier, Holmes et al. 2015
#   return(TL)
# }
# 
# # total 2 precaudal function ----
# # Function to convert total to precaudal length for various species
# # Need FAO species code and total length measure
# total_2_precaudal <- function(spp, TL) {
#   PC <- NA
#   if (spp == 'AML') { PC <- (TL - 17.88) / 1.17} # Carcharhinus amblyrhynchos, Bradley et al. 2017
#   if (spp == 'BRO') { PC <- (TL + 35.924) / 1.364 }# Carcharhinus brachyurus, Drew et al. 2016
#   if (spp == 'BSH') { PC <- 0.78 * TL - 3.75 } # Prionace glauca,  Fujinami et al. 2017
#   if (spp == 'BTH') { PC <- (TL - 0.7085) / 1.8232 }  # Alopias superciliosus, Calle-Moran et al. 2023
#   if (spp == 'CCE') { PC <- 0.81 * TL - 0.916 } # Carcharhinus leucas, Cliff & Dufley, 1991
#   if (spp == 'CCG') { PC <- 0.75865 * TL - 1.17034 }  # Carcharhinus galapagensis, Wetherbee et al. 1996
#   if (spp == 'CCP') { PC <- (TL - 4.566) / 1.316 } # Carcharhinus plumbeus, Geraghty et al. 2014
#   if (spp == 'DUS') { PC <- (TL - 8.021) / 1.305 } # Carcharhinus obscurus, Geraghty et al. 2014
#   if (spp == 'FAL') { PC <- (TL - 3.64)/1.31 } # Joung et al. 2008 PC <- (TL - 3.64) / 1.31 } # Carcharhinus falciformis, Oshitani et al. 2003
#   if (spp == 'LMA') { PC <- 0.84 * TL - 2.13 }  # Isurus paucus, Semba et al. 2009
#   if (spp == 'LMD') { PC <- (TL - 15.19) / 1.15 } # Lamna ditropis, Goldman & Musick, 2006
#   if (spp == 'OCS') { PC <- 0.755 * TL - 6.019 }  # Carcharhinus longimanus, Joung et al. 2016
#   if (spp == 'PTH') { PC <- (TL - 2.34) / 1.91 }  # Alopias pelagicus, Liu et al. 1999
#   if (spp == 'RHN') { PC <- (TL - 0.262) / 1.148 } # Rhincodon typus, Hsu et al. 2012
#   if (spp == 'SMA') { PC <- 0.84 * TL - 2.13 }  # Isurus oxyrinchus, Semba et al. 2009
#   if(spp == 'SPZ') { PC <- (TL - 4.3873) / 1.3531} # Sphyrna zygaena, Lopez Martinez et al. 2020
#   if (spp == 'TIG') { PC <- (TL - 34.321) / 1.159 } # Galeocerdo cuvier, Holmes et al. 2015
#   return(PC)
# }
# 
# # fork 2 precaudal function ----
# # Function to convert fork to precaudal length for various species
# # Need FAO species code and fork length measure
# fork_2_precaudal <- function(spp, UF) {
#   PC <- NA
#   if (spp == 'AML') { PC <- (UF - 9.78) } # Carcharhinus amblyrhynchos, Bradley et al. 2017
#   if (spp == 'BSH') { PC <-  (0.92 * UF) - 0.22 } # Prionace glauca,  Fujinami et al. 2017
#   if (spp == 'CCE') { PC <- (UF - 2.779) / 1.08 }  # Carcharhinus leucas, Cliff & Dufley, 1991
#   if (spp == 'FAL') { PC <- (UF - 1.1) / 1.09 }# Carcharhinus falciformis, Joung et al 2008
#   if (spp == 'LMA') { PC <- (0.91 * UF) - 0.95 }  # Isurus paucus, Semba et al. 2009
#   if (spp == 'LMD') { PC <- (UF - 6.91) / 1.08 } # Lamna ditropis, Goldman & Musick, 2006
#   if (spp == 'POR') { PC <- (UF - 2.1096) / 1.1068 } # Lamna nasus, Francis 2013
#   return(PC)
# }
# 
# # precaudal 2 fork function ----
# # Function to convert precaudal to fork length for various species
# # Need FAO species code and precaudal length measure
# precaudal_2_fork <- function(spp, PC) {
#   UF <- NA
#   if (spp == 'AML') { PC <- (UF - 9.78) } # Carcharhinus amblyrhynchos, Bradley et al. 2017
#   if (spp == 'BSH') { UF <- (PC + 0.22) / 0.92 } # Prionace glauca, Fujinami et al. 2017
#   if (spp == 'CCE') { UF <- 1.08 * PC + 2.779 } # Carcharhinus leucas, Cliff & Dufley, 1991
#   if (spp == 'FAL') { UF <- 1.09 * PC + 1.1}  # Carcharhinus falciformis, Joung et al 2008
#   if (spp == 'LMA') { UF <- (PC + 0.95) / 0.91 }  # Isurus paucus, Semba et al. 2009
#   if (spp == 'LMD') { UF <- 1.08 * PC + 6.91 }   # Lamna ditropis, Goldman & Musick, 2006
#   if (spp == 'POR') { UF <- 1.1068 * PC + 2.1096 } # Lamna nasus, Francis 2013
#   return(UF)
# }
# 
# # Fishery functions -----
# 
# # Selectivity function ----
# # Calculate the selectivity of fishery inc fixed, logistic, knifeedge or dome
# # Dome doesn't work currently - TBC
# selex <- function(lens, type = 'fixed', fix = 1, L1 = NA, L2 = NA, L3 = NA, m1 = NA, m2 = NA, s1 = NA, s2 = NA) {
#   if (type == 'fixed') {
#     sel <- ifelse(is.na(fix), rep(1, length(lens)), rep(fix, length(lens)))
#   } else if (type == 'knifeedge') {
#     sel <- ifelse(lens < L1, 0, 1)
#   } else if (type == 'logistic') {
#     sel <- 1 / (1 + exp(-log(19.0) * (lens - L1) / (L2 - L1)))
#     } else if (type == 'boxed') {
#     sel <- ifelse(lens < L1, 0, ifelse(lens > L2, 0, 1)) 
#   } else if (type == 'dome') {
#     if (any(is.na(c(L1, L2, L3, m1, m2, s1, s2)))) {
#       stop("Parameters L1, L2, L3, m1, m2, s1, and s2 must be provided for 'dome' type.")
#     }
#     asc <- 1 / (1 + exp(-s1)) + (1 - 1 / (1 + exp(-s1))) * (exp((-(lens - L2)^2) / m1) - exp((-(lens[1] - L2)^2) / m1)) / (1 - exp((-(lens[1] - L2)^2) / m1))
#     dsc <- 1 + (1 / (1 + exp(-s2)) - 1) * ((exp((-(lens - L3)^2) / m2) - 1) / (exp((-(lens[length(lens)] - L3)^2) / m2 - 1)))
#     J1 <- 1 / (1 + exp(-(20 * (lens - L2) / (1 + abs(lens - L2)))))
#     J2 <- 1 / (1 + exp(-20 * ((lens - L3) / (1 + abs(lens - L3)))))
#     sel <- asc * (1 - J1) + J1 * (1 - J2 + dsc * J2)
#   } else {
#     stop("Invalid type argument. Supported types are 'fixed', 'knifeedge', 'logistic', and 'dome'.")
#   }
#   return(sel)
# }
# 
# 
# easi_dome_sel <- function(lens, prop_scaled, 
#                           L1_start = NULL, L2_start = NULL, 
#                           s1_start = 10, s2_start = 10) {
#   library(zoo)
#   
#   df <- data.frame(lens = lens, prop_scaled = rollmean(prop_scaled, k = 3, fill = 0))
# 
#   
#   if (is.null(L1_start)) {
#     L1_start <- min(df$lens[df$prop_scaled >= 0.3], na.rm = TRUE)
#   }
#   if (is.null(L2_start)) {
#     L2_start <- max(df$lens[df$prop_scaled >= 0.3], na.rm = TRUE)
#   }
#   
#   # Dome-shaped selectivity function
#   dome_selectivity <- function(lens, L1, L2, s1, s2) {
#     asc <- 1 / (1 + exp(-(lens - L1) / s1))
#     dsc <- 1 - 1 / (1 + exp(-(lens - L2) / s2))
#     sel <- asc * dsc
#     sel / max(sel, na.rm = TRUE)
#   }
#   
#   # Objective function
#   obj_fun <- function(par, data) {
#     sel <- dome_selectivity(data$lens, par[1], par[2], par[3], par[4])
#     sum((data$prop_scaled - sel)^2)
#   }
#   
#   # Fit using optim
#   mod <- optim(
#     par = c(L1_start, L2_start, s1_start, s2_start),
#     fn = obj_fun,
#     data = df,
#     method = "L-BFGS-B",
#     lower = c(min(df$lens), min(df$lens) + 10, 1, 1),
#     upper = c(max(df$lens) - 10, max(df$lens), 100, 100),
#     control = list(maxit = 5000)
#   )
#   
#   # Return scaled selectivity curve
#   params <- setNames(mod$par, c("L1", "L2", "s1", "s2"))
#   sels <- dome_selectivity(df$lens, params["L1"], params["L2"], params["s1"], params["s2"])
#   
#   return(as.vector(sels))
# }
# 
# # Dome selectivity function
# dome_selectivity <- function(lens, L1, L2, s1, s2) {
#   sel <- rep(0, length(lens))
#   # Ascending limb
#   sel[lens <= L1] <- exp(-((lens[lens <= L1] - L1)^2) / (2 * s1^2))
#   # Descending limb
#   sel[lens > L1 & lens <= L2] <- 1
#   sel[lens > L2] <- exp(-((lens[lens > L2] - L2)^2) / (2 * s2^2))
#   sel <- sel / max(sel)  # Normalize to max = 1
#   return(sel)
# }
# 
# # easi_gam_sel <- function(lens, prop_scaled) {
# #   df <- data.frame(lens = lens, prop_scaled = prop_scaled)
# #   
# #   # Fit a GAM with cyclic cubic spline (or adaptive spline)
# #   mod <- tryCatch({
# #     mgcv::gam(prop_scaled ~ s(lens, bs = 'cs', k=10), data = df,
# #               family = gaussian(link = "identity"))
# #   }, error = function(e) {
# #     warning("GAM fit failed: ", e$message)
# #     return(NULL)
# #   })
# #   
# #   # If model failed, return NA vector
# #   if (is.null(mod)) {
# #     return(rep(NA, length(lens)))
# #   }
# #   
# #   # Predict from model
# #   pred <- predict(mod, newdata = df)
# #   
# #   # Clean predictions: clip and normalize
# #   pred <- pmax(pred, 0)                     # no negative selectivity
# #   pred <- pred / max(pred, na.rm = TRUE)    # scale to max 1
# #   
# #   return(as.vector(pred))
# # }
# 
# 
# 
# # easi_gam_sel2 <- function(lens, prop_scaled) {
# #   df <- data.frame(lens = lens, prop_scaled = rollmean(prop_scaled, k=3, fill = 0))
# #   
# #   best_aic <- Inf
# #   best_pred <- rep(NA, length(lens))
# #   
# #   for (k_val in 2:10) {
# #     mod <- tryCatch({
# #       mgcv::gam(prop_scaled ~ s(lens, bs = 'cs', k = k_val), data = df,
# #                 family = gaussian(link = "identity"), method = "REML", select = TRUE)
# #       #mgcv::gam(prop_scaled ~ s(lens, bs = "ts", k = k_val), data = df, method = "REML")
# #     }, error = function(e) {
# #       warning(sprintf("GAM fit failed for k=%d: %s", k_val, e$message))
# #       return(NULL)
# #     })
# #     
# #     if (!is.null(mod)) {
# #       mod_aic <- AIC(mod)
# #       if (mod_aic < best_aic) {
# #         best_aic <- mod_aic
# #         pred <- predict(mod, newdata = df)
# #         # Clean predictions: clip and normalize
# #         pred <- pmax(pred, 0)                     # no negative selectivity
# #         pred <- pred / max(pred, na.rm = TRUE)    # scale to max 1
# #         best_pred <- as.vector(pred)
# #       }
# #     }
# #   }
# #   
# #   return(best_pred)
# # }
# 
# easi_gam_sel4 <- function(lens, prop_scaled, thresh = 0.02) {
#   # Pre-smooth input with 5-point rolling mean to reduce noise
#   df <- data.frame(
#     lens = lens,
#     prop_scaled = zoo::rollmean(prop_scaled, k = 3, fill = 0)
#   )
#   
#   best_aic <- Inf
#   best_pred <- rep(NA, length(lens))
#   
#   for (k_val in 4:7) {
#     mod <- tryCatch({
#       mgcv::gam(
#         prop_scaled ~ s(lens, bs = 'cs', k = k_val),
#         data = df,
#         family = gaussian(link = "identity"),
#         method = "REML",
#         select = TRUE
#       )
#     }, error = function(e) {
#       warning(sprintf("GAM fit failed for k = %d: %s", k_val, e$message))
#       return(NULL)
#     })
#     
#     if (!is.null(mod)) {
#       mod_aic <- AIC(mod)
#       
#       if (mod_aic < best_aic) {
#         best_aic <- mod_aic
#         
#         # Predict and clean
#         pred <- predict(mod, newdata = df)
#         pred <- pmax(pred, 0)
#         pred <- pred / max(pred, na.rm = TRUE)
#         
#         # --- Enforce dome shape: monotonic taper from peak ---
#         peak_index <- which.max(pred)
#         
#         # Left side: ensure non-increasing toward the left
#         if (peak_index > 1) {
#           for (i in (peak_index - 1):1) {
#             if (pred[i] > pred[i + 1]) {
#               pred[i] <- pred[i + 1]
#             }
#           }
#         }
#         
#         # Right side: ensure non-increasing toward the right
#         if (peak_index < length(pred)) {
#           for (i in (peak_index + 1):length(pred)) {
#             if (pred[i] > pred[i - 1]) {
#               pred[i] <- pred[i - 1]
#             }
#           }
#         }
#         
#         best_pred <- as.vector(pred)
#       }
#     }
#   }
#   
#   return(best_pred)
# }
# 
# 
# 
# # easi_dome_sel <- function(lens, prop_scaled, s1_start = 10, s2_start = 10) {
# #   
# #   df <- data.frame(lens = lens, prop_scaled = prop_scaled)
# #   # Estimate L1 and L2 from where prop_scaled is near 0.4
# #   L1_start <- df$lens[which.max(df$prop_scaled[df$prop_scaled <= 0.4])]
# #   L2_start <- df$lens[which.max(df$prop_scaled[df$prop_scaled <= 0.6])]
# #   
# # # Define your dome-shaped selectivity function
# # dome_selectivity <- function(lens, L1, L2, s1, s2) {
# #   asc <- 1 / (1 + exp(-(lens - L1) / s1))
# #   dsc <- 1 - 1 / (1 + exp(-(lens - L2) / s2))
# #   sel <- asc * dsc
# #   sel / max(sel, na.rm = TRUE)
# # }
# # 
# # # Objective function for minimization
# # obj_fun <- function(par, data) {
# #   sel <- dome_selectivity(data$lens, par[1], par[2], par[3], par[4])
# #   sum((data$prop_scaled - sel)^2)
# # }
# # 
# # mod <- optim(
# #   par = c(L1_start, L2_start, s1_start, s2_start),
# #   fn = obj_fun,
# #   data = df,
# #   method = "L-BFGS-B",
# #   lower = c(min(df$lens), min(df$lens) + 10, 1, 1),
# #   upper = c(max(df$lens) - 10, max(df$lens), 100, 100),
# #   control = list(maxit = 5000)
# # )
# # 
# # params <- setNames(mod$par, c("L1", "L2", "s1", "s2"))
# # 
# # sels <- dome_selectivity(lens, params["L1"], params["L2"], params["s1"], params["s2"])
# # return(as.vector(sels))
# # 
# # }
# 
# # Easi functions ----
# 
# # encounterability function ----
# # this calculates species encounterability based on the vertical distribution of spp vs fishery
# encount <- function(min_spp, max_spp, min_gr, max_gr) {
#   if (max_spp <= min_gr | min_spp >= max_gr) { 
#     enc <- 0
#   } else if (max_spp <= max_gr & min_spp >= min_gr) { 
#     enc <- 1
#   } else { 
#     enc <- (max(min(max_gr, max_spp) - max(min_gr, min_spp), 0)) / (max_spp - min_spp)
#   }
#   return(enc)
# }
# 
# 
# # cumulative product function ----
# # This is to match excel product function
# # Handles 0s and assigns values at right row index
# cumprod2 <- function(x) {
#   result <- rep(0, length(x))
#   nz_indices <- which(x != 0)
#   result[nz_indices] <- cumprod(x[nz_indices])
#   return(result)
# }
# 
# # Natural mortality (M) functions ----
# 
# # Griffiths EASIfish M function ----
# # This function recreates the Griffiths et al 2019 M estimator flow diagram
# # This uses a range of M functions (defined below) to estimate M based on data available.
# # It prioritises a direct estimate M, then tmax based estimators, Linf and others.
# griffiths_M_estimator <- function(M = NA, tmax = NA, Linf = NA, K = NA, sst = NA, Lmax = NA) {
#   
# # Calculate all possible M estimates
#   M_hoenig83 <- hoenig83(tmax)
#   M_then_Tmax <- then_Tmax(tmax)
#   M_then_Linf <- then_Linf(K = K, Linf = Linf)
#   M_pauly_kt <- pauly_kt(sst = sst, K = K)
#   M_jensen <- jensen(K = K)
#   M_pauly_lt <- pauly_lt(sst = sst, Linf = Linf, Lmax = Lmax)
#   #Linf <- Lmax_to_Linf(Lmax)
#   M1 <- mean(c(M_hoenig83, M_then_Tmax))
#   M2 <- mean(c(M_then_Linf, M_pauly_kt))
#   M3 <- mean(c(M1, M2))
#   M4 <- mean(c(M_pauly_kt, M_jensen))
#   M5 <- mean(c(M4, M_pauly_lt))
#   M_pauly_lt2 <- pauly_lt(Linf = NA, sst = sst, Lmax = Lmax)
#   
#   # Now apply flow diagram logic based on which estimates are available
#   if (!is.na(M)) {
#     out <- M
#   } else if (!is.na(M3)) {
#     out <- M3
#   } else if (!is.na(M1) & is.na(M2)) {
#     out <- M1
#   } else if (is.na(M1) & !is.na(M2)) {
#     out <- M2
#   } else if (is.na(M1) & is.na(M2) & !is.na(M5)) {
#     out <- M5
#   } else if (is.na(M1) & is.na(M2) & is.na(M5) & !is.na(Linf)) {
#     out <- M_pauly_lt
#   } else if (is.na(M1) & is.na(M2) & is.na(M5) & is.na(Linf) & !is.na(K)) {
#     out <- M4
#   } else if (is.na(M1) & is.na(M2) & is.na(M5) & is.na(K) & is.na(Linf)) {
#     out <- M_pauly_lt2
#   }
#   return(out)
# }
# 
# # griffiths_M_estimator(M=NA, tmax = 20, K = 0.1, Linf = 200, sst = 27, Lmax = 210)
# 
# # easiM_estimator function ----
# # Updated M estimator for SPC sharks easifish 2025
# # Function to extract a mean M with lognormal lower and upper estimates.
# # Estimates M from a vector of function names, takes the mean and estimates CI around them based on CV value
# # For this study, the mean of 3 tmax estimators were used to derive M values
# 
# # Choose from: 
# #hoenig83(tmax = tmax), 
# #hoenig_fish(tmax = tmax),
# #jensen(K = K),
# #pauly_sst(sst = sst, Linf = Linf, K = K),
# #pauly_kt(sst = sst, K = K), 
# #durueil(tmax = tmax),
# #durueil_shark(tmax = tmax), 
# #durueil_Linf(Linf = Linf, K = K, L0 = L0, X = 0.99),
# #then_Tmax(tmax = tmax), 
# #then_Linf(K = K, Linf = Linf),
# #frisk1(K), 
# #frisk2(tmat), 
# #Liu_Linf(K = K, Linf = Linf, tmax = NA, method = 'Linf'),
# #Liu_tmax(K = K, Linf = Linf, tmax = tmax, method = 'tmax')
# # args = K, Linf, tmax, temp, L0, CV
# 
# easiM_estimate <- function(funcs, ...) {
#   args <- list(...)
#   
#   # Get individual M estimates
#   M_vals <- sapply(funcs, function(f_name) {
#     f <- tryCatch(match.fun(f_name), error = function(e) return(NA))
#     if (is.function(f)) {
#       f_args <- names(formals(f))
#       usable_args <- args[names(args) %in% f_args]
#       result <- tryCatch(do.call(f, usable_args), error = function(e) NA)
#       return(as.numeric(result))  # Ensure the result is numeric
#     } else {
#       return(NA)
#     }
#   })
#   
#   mean_M <- as.vector(mean(M_vals, na.rm = TRUE))
#   
#   # Return as a data frame
#   return(mean_M)
# }
# 
# easiM_estimator <- function(funcs, CV = 0.2, ...) {
#   args <- list(...)
#   
#   # Get individual M estimates
#   M_vals <- sapply(funcs, function(f_name) {
#     f <- tryCatch(match.fun(f_name), error = function(e) return(NA))
#     if (is.function(f)) {
#       f_args <- names(formals(f))
#       usable_args <- args[names(args) %in% f_args]
#       tryCatch(do.call(f, usable_args), error = function(e) NA)
#     } else {
#       NA
#     }
#   })
#   
#   mean_M <- mean(M_vals, na.rm = TRUE)
#   
#   # Handle edge case
#   if (is.na(mean_M) || mean_M <= 0) {
#     return(data.frame(mean_M = NA, lower = NA, upper = NA))
#   }
#   
#   # Lognormal distribution parameters from mean and CV
#   sdlog <- sqrt(log(1 + CV^2))
#   meanlog <- log(mean_M) - (sdlog^2) / 2
#   
#   # 95% confidence interval from lognormal
#   lower <- qlnorm(0.025, meanlog, sdlog)
#   upper <- qlnorm(0.975, meanlog, sdlog)
#   
#   # Return as a data frame
#   return(data.frame(M_est = mean_M, M_lower = lower, M_upper = upper))
# }
# 
# # easi_M_estimator(funcs = c("hoenig83", "Liu_tmax", "Liu_Linf", 
# #                            "jensen", "pauly_sst","pauly_kt", "durueil",
# #                            "durueil_shark", "hoenig_fish", "durueil_Linf", "frisk1",
# #                            "Then_tmax", "Then_Linf"), CV = 0.1, tmax = 20, K = 0.1, L0 = 60)
# 
# # easiM_estimator_all function ----
# # Same function as above, but returns a df of all M values rather than a mean with bounds, plus optional plot
# 
# easiM_estimator_all <- function(funcs, plot = FALSE, ...) {
#   args <- list(...)
#   
#   results <- lapply(funcs, function(f_name) {
#     f <- tryCatch(match.fun(f_name), error = function(e) return(NULL))
#     if (!is.null(f)) {
#       f_args <- names(formals(f))
#       usable_args <- args[names(args) %in% f_args]
#       
#       M_val <- tryCatch(do.call(f, usable_args), error = function(e) NA)
#       
#       if (!is.na(M_val)) {
#         return(data.frame(func_name = f_name, M = M_val))
#       }
#     }
#     return(NULL)
#   })
#   
#   # Combine valid rows into a data frame
#   M_df <- do.call(rbind, results)
#   
#   # Add row with mean M
#   if (!is.null(M_df) && nrow(M_df) > 0) {
#     mean_M <- mean(M_df$M, na.rm = TRUE)
#     M_df <- rbind(M_df, data.frame(func_name = "mean", M = mean_M))
#   } else {
#     M_df <- data.frame(Estimator = "mean", M_est = NA)
#   }
#   # Optional plot
#   if (plot) {
#     p <- ggplot(M_df[M_df$func_name != "mean", ], aes(x = func_name, y = M)) +
#       geom_point(color = "blue", size = 3) +
#       geom_hline(yintercept = M_df$M[M_df$func_name == "mean"], linetype = "dashed", color = "red") +
#       labs(#title = "Natural Mortality Estimates by Function",
#            x = "Estimator",
#            y = "Natural mortality") +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       ylim(0,NA)
#     
#     print(p)
#   }
#   return(M_df)
# 
# }
# 
# # easi_M_estimator_all(funcs = c("hoenig83", "Liu_tmax", "Liu_Linf", 
# #           "jensen", "pauly_sst","pauly_kt", "durueil",
# #            "durueil_shark", "hoenig_fish", "durueil_Linf", "frisk1",
# #            "Then_tmax", "Then_Linf"), CV = 0.1, tmax = 20, K = 0.1, L0 = 60, plot = T)  
# 
# # Lmax to Linf function ----
# # Convert an Lmax estimate to Linf from Froese and Binolahn 2000
# Lmax_to_Linf <- function(Lmax) {
#   Linf <- 10 ^ (0.044 + 0.9841 * log(Lmax, base = 10))
#   return(Linf)
# }
# 
# # Hoenig 1983 M function ----
# # Estimate M function from Hoenig 1983
# hoenig83 <- function(tmax) {
#   M <- 4.3/tmax
#   return(M)
# }
# 
# # Hoenig M function for fishes ----
# hoenig_fish <- function(tmax) {
#   M <- exp(1.46 - 1.01 * log(tmax))
#   return(M)
# }
# 
# # Pauly M function with Linf and SST ----
# pauly_lt <- function(Linf, sst, Lmax) {
#   if(!is.na(Linf)) {M <- 10 ^ (0.556 - 0.718 * log(Linf)) + (0.02 * sst)}
#   if(is.na(Linf)) {Linf <- 10 ^ (0.044 + 0.9841 * log(Lmax, base = 10))
#   M <- 10 ^ (0.556 - 0.718 * log(Linf)) + (0.02 * sst)}
#   return(M)
# }
# 
# # Pauly M function with Linf, K and sst ----
# pauly_sst <- function(sst, Linf, K) {
#   if (sst < 0.001 || is.na(Linf) || is.na(K) || Linf == "" || K == "") {
#     return("")
#   } else {
#     M <- exp(-0.0066 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(sst))
#     return(M)
#   }
# }
# 
# # Pauly M function with K and sst ----
# pauly_kt <- function(sst, K) {
#   if (is.na(sst) || is.na(K) || sst == "" || K == "") {
#     return("")
#   } else {
#     M <- K * exp(-0.22 + 0.3 * log(sst))
#     return(M)
#   }
# }
# 
# # jensen M function ----
# jensen <- function(K) {
#   if(is.na(K)) { return("")}
#   else{
#     M <- 1.6 * K
#     return(M)
#   }
# }
# 
# # Durueil et al. 2021 M function for all species ----
# durueil <- function(tmax) {
#   if(is.na(tmax) | !is.numeric(tmax)) { return("Error: Check tmax")}
#   M <- exp(1.551 - 1.061 * log(tmax))
#   return(M)
# }
# 
# # Durueil et al. 2021 M function for sharks ----
# durueil_shark <- function(tmax) {
#   if(is.na(tmax) | !is.numeric(tmax)) { return("Error: Check tmax")}
#   M <- exp(1.583 - 1.087 * log(tmax))
#   return(M)
# }
# 
# # Durueil et al. 2021 M function with Linf ----
# durueil_Linf <- function(Linf, L0, K, X = 0.99) {
#   get_Tmax <- function(Linf, K, L0, X) {
#     tmax <- (1/K) * log((Linf - L0)/((1 - X) * L0))
#     return(tmax)
#   }
#   M <- exp(1.583 - 1.087 * log(get_Tmax(Linf, K, L0, X)))
#   return(M)
# }
# 
# # Then et al. 2015 M function with T max also known as Hoenig_nls ----
# then_Tmax <- function(tmax) {
#   M <- 4.899 * tmax ^ -0.916
#   return(M)
# }
# 
# # Then et al. 2015 M function with Linf also known as Pauly_nls ----
# then_Linf <- function(K, Linf) {
#   M <- 4.118 * (K ^ 0.73) * (Linf ^ -0.33)
#   return(M)
# }
# 
# # Frisk et al 2001 M function for sharks with K ----
# frisk1 <- function(K) {
#   M <- 0.42 * log(K) - 0.83
#   M <- exp(M)
#   return(M)
# }
# 
# # Frisk et al 2001 M function for sharks with age at maturity ----
# frisk2 <- function(tmat) {
#   M <- 1 / (0.44 * tmat + 1.87)
#   return(M)
# }
# 
# # Liu_tmax function ----
# # Liu et al 2021 natural mortality function using BRT
# # Note, you need metadata R data file.
# # get from here - https://github.com/ChanjuanLiu92/M-estimate
# 
# Liu_tmax <- function(K = NA, Linf = NA, tmax = NA, Liu_method = 'tmax'){
#   library(caret)
#   library(gbm)
# 
#   if(Liu_method == 'K') {Liu_method <- 'BRT1'}
#   if(Liu_method == 'Linf') {Liu_method <- 'BRT2'}
#   if(Liu_method == 'tmax') {Liu_method <- 'BRT3'}
#   
#   # Deifne data to predict to
#   mydata <- data.frame(K = K, Linf = Linf, tmax = tmax, Class = 1)
#   
#   # Create data to model to
#   load('metadata324.Rdata')
#   M=metadata324$M
#   K=metadata324$K
#   Linf=metadata324$Linf
#   tmax=metadata324$tmax
#   Class=as.factor(metadata324$Class)
#   
#   data1=na.omit(data.frame(M,tmax,Class))
#   data2=na.omit(data.frame(M,K,Linf,Class))
#   data3=na.omit(data.frame(M,K,Linf,tmax,Class))
#   
#   if (Liu_method=="BRT1"){
#     model1=gbm(M~.,data=data1,distribution="gaussian",
#                n.trees=400,
#                interaction.depth=5,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model1,mydata,n.trees=400))[1,]
#   }
#   
#   else if (Liu_method=="BRT2"){
#     model2=gbm(M~.,data=data2,distribution="gaussian",
#                n.trees=300,
#                interaction.depth=8,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model2,mydata,n.trees=300))[1,]
#   }
#   
#   else if (Liu_method=="BRT3"){
#     model3=gbm(M~.,data=data3,distribution="gaussian",
#                n.trees=380,
#                interaction.depth=8,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model3,mydata,n.trees=380))[1,]
#   }
#   return(out)
# }
# 
# # Liu_Linf function ----
# # Liu et al 2021 natural mortality function using BRT
# # Note, you need metadata R data file.
# # get from here - https://github.com/ChanjuanLiu92/M-estimate
# 
# Liu_Linf <- function(K = NA, Linf = NA, tmax = NA, Liu_method = 'Linf'){
#   library(caret)
#   library(gbm)
# 
#   if(Liu_method == 'K') {Liu_method <- 'BRT1'}
#   if(Liu_method == 'Linf') {Liu_method <- 'BRT2'}
#   if(Liu_method == 'tmax') {Liu_method <- 'BRT3'}
#   
#   # Define data to predict to
#   mydata <- data.frame(K = K, Linf = Linf, tmax = tmax, Class = 1)
#   
#   # Create data to model to
#   load('metadata324.Rdata')
#   M=metadata324$M
#   K=metadata324$K
#   Linf=metadata324$Linf
#   tmax=metadata324$tmax
#   Class=as.factor(metadata324$Class)
#   
#   data1=na.omit(data.frame(M,tmax,Class))
#   data2=na.omit(data.frame(M,K,Linf,Class))
#   data3=na.omit(data.frame(M,K,Linf,tmax,Class))
#   
#   if (Liu_method=="BRT1"){
#     model1=gbm(M~.,data=data1,distribution="gaussian",
#                n.trees=400,
#                interaction.depth=5,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model1,mydata,n.trees=400))[1,]
#   }
#   
#   else if (Liu_method=="BRT2"){
#     model2=gbm(M~.,data=data2,distribution="gaussian",
#                n.trees=300,
#                interaction.depth=8,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model2,mydata,n.trees=300))[1,]
#   }
#   
#   else if (Liu_method=="BRT3"){
#     model3=gbm(M~.,data=data3,distribution="gaussian",
#                n.trees=380,
#                interaction.depth=8,
#                shrinkage=0.01,
#                n.minobsinnode = 1)
#     out <- data.frame(predict(model3,mydata,n.trees=380))[1,]
#   }
#   return(out)
# }
# 
# # Mage_ChenWat function ----
# # M at age by Chen and Watanabe 1989
# # Chen, S. and S. Watanabe. 1989. Age Dependence of Natural Mortality Coefficient in Fish Population Dynamics. Nippn Suisan Gakkaishi 55(2): 205-208.
# # These are pulled from Cope et al 2021 R shiny app code
# Mage_ChenWat <- function(Age_in,K,t0)
# {
#   if(anyNA(c(Age_in,K,t0))){M.out<-NA}
#   else
#   {
#     a <- Age_in
#     tM <- -1/K * (log(abs(1 - exp(K * t0))) + t0)
#     a0 <- 1 -exp(-K * (tM - t0))
#     a1 <- K * exp(-K * (tM - t0))
#     a2 <- -0.5 * K^2 * exp(-K * ( tM - t0))
#     if(a <= tM) {M.out <- K / (1 - exp(-K * (a - t0)))}
#     if(a > tM) {M.out <- K / (a0 + a1 * (a - tM) + a2 * (a - tM)^2)}
#   }
#   return(M.out)
# }
# 
# # Mage_Gislason2 function ----
# # M at age Gislason 2010
# # Gislason, H., N. Daan, J. C. Rice, and J. G. Pope. 2010. Size, growth, temperature and the natural mortality of marine fish. Fish and Fisheries 11: 149-158.
# 
# # Original from Cope that creates lengths then runs M at age function
# # using M.empirical() from fishmethods
# # Mage_Gislason <- function(Amax,Linf,k,t0) {
# #   Lts<-Linf*(1-exp(-k*(c(1:Amax)-t0)))
# #   Gis_Ms_a<-mapply(function(x) M.empirical(Linf=Linf,Kl=k,Bl=Lts[x],method=9)[1],x=1:length(Lts),SIMPLIFY=TRUE)
# #   return(Gis_Ms_a)
# # }
# 
# # Updated function for use within easifish code/data structure
# Mage_Gislason2 <-function(len,Amax,Linf,K,t0) {
#   Gis_Ms_a<- M.empirical(Linf=Linf,Kl=K,Bl=len,method=9)[1] 
#   return(Gis_Ms_a)
# }
# 
# # Mage_charnov2 function ----
# # M at age Charnov 2013
# # Charnov, E.L., Gislason, H., Pope, J.G., 2013. Evolutionary assembly rules for fish life histories. Fish and Fisheries 14, 213-224. https://doi.org/10.1111/j.1467-2979.2012.00467.x
# 
# # Original from Cope et al 2021 using M.empirical() from fishmethods
# # Mage_charnov <- function(Amax,Linf,K,t0) {
# #   Lts<-Linf*(1-exp(-K*(c(1:Amax)-t0)))
# #   Charnov_Ms_a<-mapply(function(x) M.empirical(Linf=Linf,Kl=K,L=Lts[x],method=13)[1],x=1:length(Lts),SIMPLIFY=TRUE)
# #   return(Charnov_Ms_a)
# # }
# 
# # Updated function to fit easifish data structure
# Mage_charnov2 <- function(len, Amax,Linf,K,t0) {
#   Charnov_Ms_a <- M.empirical(Linf=Linf,Kl=K,L=len,method=13)[1]
#   return(Charnov_Ms_a)
# }
# 
# #Mage_charnov2(lens = 200, Amax = 30,Linf = 211, K = 0.09, t0 = -6.1)
# #map_dbl(seq(10, 300, by = 10), ~ Mage_charnov2(.x, Amax = 30, K = 0.06, t0 = -6, Linf = 211))
# 
# # x <- M_estimator(tmax = 22, K = 0.129, Linf = 251.9, sst = 29)
# # x |>
# #   dplyr::filter(Function != 'mean') |>
# #   ggplot() +
# #   aes(x = factor(0), y = Estimate) +
# #   geom_boxplot()
# 
# ####################################################
# # Not currently used
# 
# # Then 2015 updated Hoenig 1983 function ----
# # Updated Hoenig 1983 function from Then et al. 2015
# # hoenig_nls <- function(tmax) {
# #   M <- 4.899 * tmax ^ -0.916
# #   return(M)
# # }
# 
# # Then 2015 updated Pauly Linf function ----
# # Updated Pauly function from Then et al. 2015
# # pauly_nls <- function(K, Linf) {
# #   M <- (4.118 * K^0.73) * (Linf^-0.33)
# #   return(M)
# # }
# 
# ############################################
# # slope function
# # slope - calculates the slope between x1/x2 and y1/y2 from lm. Used for ypr calculations
# #slope <- function(x, y) cov(x, y) / var(x)
# #intercept <- function(x, y) mean(y) - slope(x, y) * mean(x)
# 
# ###########################################
# 
# # SVGBF function
# # Von Bertalanffy growth function estimating age from length
# # SVGBF <- function(len, Linf, K, t0) {
# #   age <- log(1 - (len/Linf)) / -K + t0
# #   return(age)
# # }
# 
# # SVGBF function
# # Von Bertalanffy growth function estimating age from length
# # GVGBF <- function(len, Linf, K, t0, D) {
# #   age <- log(1 - (len/Linf)^D)  / -K + t0
# #   return(age)
# # }
# 
# # VBL0 function
# # Von Bertalanffy growth function estimating age from length zero for sharks
# # VBL0 <- function(len, Linf, K, L0) {
# #   age <- 1 / K * log((Linf - L0) / (Linf - len))
# #   return(age)
# # }
# 
# # Gompertz growth function
# # Estimate growth using gompertz function
# # Gompertz <- function(len, Linf, K, t0) {
# #   age <- -log(-log(len / Linf) / (1/K)) / K + t0
# #   return(age)
# # }
# 
# # Logistic growth function
# # Estimate growth using logistic growth equation
# # Logistic <- function(len, Linf, K, t0) {
# #   age <- t0 + (-log((Linf - len)/len) / K)
# #   return(age)
# # }
# 
# # Estimate growth using richards growth equation
# # Richards <- function(len, Linf, K, t0, m) {
# #   age <- -log(((Linf/len)^m-1) / t0) / K
# #   return(age)
# # }
# 
# # Gompertz <- function(len, Linf, gi = NULL, ti = NULL) {
# #   # gi = growth at inflection point
# #   # ti = age at inflection pt
# #   # see FSA:::gompertzFuns
# #   age <- Linf * exp(-exp(-gi * (t - ti)))
# #   return(age)
# # }
# # 
# # logistic <- function(t, Linf, gninf = NULL, ti = NULL) {
# #   # ti = age at inflection pt
# #   # gninf = growth rate at -infinity?
# #   # see FSA:::logisticFuns
# #   age <- Linf/(1 + exp(-gninf * (t - ti)))
# #   return(age)
# # }
# # 
# # richards <- function(age, Linf, K, a, b) {
# #   # K = slope at max growth rate
# #   # b = length (y axis) of inflection pt
# #   # a = age (x axis) of inflection pt
# #   # s = inflection  len pt (y axis)
# #   len <- Linf * (1 - a * exp(-K * t))^b
# #   len <- Linf * (1 + (s - 1) * exp(-K(age - y))) ^ (1/(1 - s))
# #   return(len)
# # }
# 
# # schnute <- function(age, len, Linf, K, a1) {
# #   #len <- Linf + (len1 - Linf) * exp(-K(age-age1))
# #   age <- age1 - (1/K) * log((len/L0)/(len1 - Linf))
# #   return(age)
# # }
# # 
# # Schnute2 <- function(L1, L2, K, len, a, b, t1, t2) {
# #   #len <- (L1^b + (L2^b - L1^b) * ((1 - exp(-a * (age-t1)))/(1 - exp(-a * (t2-t1)))))^(1/b)
# #   age <- t1 + log(1 - ((len^b - L1^b) / (L2^b - L1^b)) * (1 - exp(-a * (t2 - t1)))) / (-a)
# # return(age)
# # }
# # 
# # x <- data.frame(len = c(1:800)) |>
# #   rowwise() |>
# #   mutate(age = Schnute2(L1=151.8, L2=466.82, t1=1, t2=50, a=0.65, b=-19.23, len = len))
# # 
# # x |> 
# #   ggplot() +
# #   aes(age, len) +
# #   geom_line() +
# #   xlim(0,50)
# # 
# # Gompertz <- function(len, Linf, L0, K) {
# #   #len <- L0 * exp(log(Linf/L0) * (1 - exp(-K * age)))
# #   age <- -1 / K * log(1 - log(len / L0) / log(Linf / L0))
# #   return(age)
# # }
# 
# ######################################################
# # Linf from Lmax
# 
# # Lmax2Linf <- function(Lmax) {
# #   Linf <- 10^(.044 + 0.9841 * log10(Lmax))
# #   return(Linf)
# # }
# # 
# # Phi2K <- function(phi, Linf) {
# #   if (is.na(phi) || is.na(Linf) || phi == "" || Linf == "") {
# #     return(NA)  # or return an appropriate value for empty or missing input
# #   } else {
# #     K <- 10^(phi - 2 * log10(Linf))
# #     return(K)
# #   }
# # }
# # 
# # TL2FL <- function(a, b, TL) {
# # FL <- a +(b * TL)
# # return(FL)
# #   }
# # 
# # Linf2Lm <- function(Linf, CI_lower = NA, CI_upper = NA, CI_level = 95) {
# #    Lm <- 10^(.9459 * log10(Linf) - 0.165)
# #    return(Lm)
# # }
# # 
# # ######################################
# # # Post capture mortality - not in excel sheet
# # # equation in griffiths et al 2023
# # # This is by length class and fishery
# # pcm_func <- function(avm, prm) {
# #   if(is.na(avm) | is.na(prm)) {return(NA)}
# #   else {pcm <- avm + ((1 - avm) * prm)
# #   return(pcm)}
# # }
# 
# 
