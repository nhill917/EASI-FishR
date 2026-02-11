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
