################################################################################
## succession_traits: functions.R
## Helper functions sourced by all analysis scripts (02–06).
##
## Sections:
##   1. PCA helpers                    (used by 02_environment_pca.R)
##   2. RF modelling helpers           (used by 03_rf_fit.R, 05_pdp.R)
##   3. SHAP helpers                   (used by 04_shap.R)
##   4. Partial dependence helpers     (used by 05_pdp.R)
##   5. Predictability helpers         (used by 06_vecv.R)
##
## Author: M. Weiss @ Maynard Lab UCL / ETH Zürich
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# 1. PCA helpers
# ══════════════════════════════════════════════════════════════════════════════

# Filter varimax-rotated scores to retain only variables with loadings above
# a given threshold. Used in 02_environment_pca.R for diagnostic display only
# (the actual predictors use full rotated scores, not filtered ones).
filter_scores <- function(scores, loadings, threshold = 0.2) {
	strong  <- abs(loadings) >= threshold
	out     <- matrix(NA_real_, nrow = nrow(scores), ncol = ncol(scores))
	for (i in seq_len(ncol(loadings))) {
		out[, i] <- ifelse(strong[, i], scores[, i], NA_real_)
	}
	rownames(out) <- rownames(scores)
	colnames(out) <- colnames(scores)
	out
}

# Compute per-variable contributions (%) to each rotated component.
# Input: varimax loadings matrix. Output: matrix of same dimensions.
varimax_contrib <- function(loadings) {
	loadings^2 / colSums(loadings^2) * 100
}

# Barplot of variable contributions to one PC, with abbreviated axis labels.
create_contrib_plot <- function(var_contrib_data, title, fill_color, abbreviations) {
	var_contrib_data %>%
		tibble::as_tibble() %>%
		dplyr::mutate(name = abbreviations[match(name, names(abbreviations))]) %>%
		ggplot2::ggplot(ggplot2::aes(x = reorder(name, -contrib), y = contrib)) +
		ggplot2::geom_bar(stat = "identity", colour = "black",
											fill = fill_color, alpha = 0.8) +
		ggplot2::labs(y = "Contributions (%)", x = NULL, title = title) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			text         = ggplot2::element_text(family = "sans", size = 8),
			axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
		)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. RF modelling helpers
# ══════════════════════════════════════════════════════════════════════════════

# Split a data frame into train/test sets by random row sampling.
# Returns a named list: list(train = ..., test = ...).
split_train_test <- function(df, ntest = NULL, test_frac = 0.2, seed = 42L) {
	set.seed(seed)
	df    <- dplyr::mutate(df, .id = dplyr::row_number())
	ntest <- if (is.null(ntest)) max(1L, floor(nrow(df) * test_frac)) else
		min(as.integer(ntest), nrow(df))
	test_ids <- sample(df$.id, size = ntest, replace = FALSE)
	list(
		train = dplyr::filter(df, !.id %in% test_ids),
		test  = dplyr::filter(df,  .id %in% test_ids)
	)
}

# Grid search over hyperparameters using OOB error on training data.
# Parallelises across grid rows with foreach %dopar%.
# Returns: list(error_plot, best_hyperparameters).
# NOTE: column names in hyper_grid use dot notation (num.trees, min.node.size)
# because expand.grid() produces dots. best_hyperparameters renames to
# underscores for compatibility with fit_rf_model().
tune_rf_model <- function(trait, data, covariates, hyper_grid, num_threads = 1L) {
	
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	hg      <- hyper_grid
	
	pred_err <- foreach::foreach(
		i           = seq_len(nrow(hg)),
		.combine    = "c",
		.packages   = "ranger"
	) %dopar% {
		row <- hg[i, ]
		mod <- ranger::ranger(
			formula    = formula, data = data,
			num.trees  = row$num.trees,
			mtry       = row$mtry,
			min.node.size = row$min.node.size,
			num.threads   = 1L,
			importance    = "none",
			respect.unordered.factors = "order",
			keep.inbag    = FALSE,
			seed          = 15L
		)
		mod$prediction.error
	}
	
	error_plot <- dplyr::bind_cols(hg, prediction_error = pred_err) %>%
		ggplot2::ggplot(
			ggplot2::aes(x = mtry,
									 y = as.factor(min.node.size),
									 fill = prediction_error)
		) +
		ggplot2::facet_wrap(~ num.trees) +
		ggplot2::geom_tile() +
		ggplot2::scale_y_discrete(breaks = c(1, 10, 20)) +
		ggplot2::scale_fill_viridis_c() +
		ggplot2::ylab("min.node.size") +
		ggplot2::ggtitle(trait) +
		ggplot2::theme(
			text              = ggplot2::element_text(family = "sans", size = 6),
			legend.key.height = grid::unit(3, "mm"),
			legend.key.width  = grid::unit(3, "mm"),
			plot.title        = ggplot2::element_text(face = "bold", size = 6)
		)
	
	best <- dplyr::bind_cols(hg, prediction_error = pred_err) %>%
		dplyr::mutate(trait = trait) %>%
		dplyr::arrange(prediction_error) %>%
		dplyr::slice(1) %>%
		# Rename to underscore convention expected by fit_rf_model()
		dplyr::rename(num_trees = num.trees, min_node_size = min.node.size)
	
	list(error_plot = error_plot, best_hyperparameters = best)
}

# Fit a single ranger RF model for one trait using pre-tuned hyperparameters.
# hyper_parameters must have columns: trait, num_trees, mtry, min_node_size.
# Returns: list(trait_mod = ranger object, performance = tibble).
fit_rf_model <- function(trait, df_train, covariates, hyper_parameters,
												 num_threads = 1L, case_weights_col = NULL) {
	
	stopifnot(
		"trait not found in training data" = trait %in% names(df_train)
	)
	
	formula <- as.formula(paste(trait, "~", paste(covariates, collapse = " + ")))
	
	row <- hyper_parameters[hyper_parameters$trait == trait, , drop = FALSE]
	if (nrow(row) != 1L)
		stop(sprintf("Hyper-grid has %d rows for trait '%s' (expected 1).",
								 nrow(row), trait))
	
	cw <- if (!is.null(case_weights_col) &&
						case_weights_col %in% names(df_train))
		df_train[[case_weights_col]] else NULL
	
	mod <- ranger::ranger(
		formula       = formula,
		data          = df_train,
		num.trees     = row$num_trees[1L],
		mtry          = row$mtry[1L],
		min.node.size = row$min_node_size[1L],
		num.threads   = num_threads,
		case.weights  = cw,
		importance    = "none",
		respect.unordered.factors = "order",
		keep.inbag    = FALSE,
		seed          = 15L
	)
	
	perf <- tibble::tibble(
		trait         = trait,
		mtry          = mod$mtry,
		num_trees     = mod$num.trees,
		min_node_size = mod$min.node.size,
		rsq           = mod$r.squared,
		pred_error    = mod$prediction.error
	)
	
	list(trait_mod = mod, performance = perf)
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. SHAP helpers
# ══════════════════════════════════════════════════════════════════════════════

# Prediction wrapper for fastshap::explain().
# ranger returns $predictions; this wrapper ensures a plain numeric vector.
predict_fn <- function(object, newdata) {
	p <- predict(object, data = as.data.frame(newdata))
	if (!is.null(p$predictions)) return(p$predictions)
	if (!is.null(p$prediction))  return(p$prediction)
	stop("Unknown predict() return structure from ranger.")
}

# Robust symmetric scaling to [-1, 1] for SHAP beeswarm colour mapping.
# Clips at ±3 SD to prevent extreme values from washing out the colour scale.
symmetric_scale <- function(x) {
	x <- as.numeric(x)
	if (!any(is.finite(x))) return(rep(0, length(x)))
	z <- as.numeric(scale(x))
	z[!is.finite(z)] <- 0
	pmax(pmin(z, 3), -3) / 3
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. Partial dependence helpers
# ══════════════════════════════════════════════════════════════════════════════

# Stratify a data frame into lower and upper quantile groups for one variable.
# Returns: list(lower = ..., upper = ...).
# Uses .data[[variable]] for safe tidy-eval inside dplyr::filter().
stratify_within <- function(data, variable, probs = c(0.25, 0.75)) {
	qs <- quantile(data[[variable]], probs = probs, na.rm = TRUE)
	list(
		lower = dplyr::filter(data, .data[[variable]] <= qs[1]),
		upper = dplyr::filter(data, .data[[variable]] >= qs[2])
	)
}

# Compute a marginal partial dependence curve for one model and one feature.
# Averages predictions over all observations while varying the focal feature
# across an evenly-spaced grid of n_grid values within its observed range.
#
# This replaces the old calculate_partial_dependence() which:
#   - called pdp::partial() redundantly then overwrote its output
#   - used foreach %dopar% internally, causing nested parallelism when called
#     from within a bootstrap foreach loop
#
# Returns: tibble with columns <feature> (grid values) and yhat (mean prediction).
compute_pdp <- function(model, data, feature, n_grid = 20L) {
	
	grid_vals <- seq(
		min(data[[feature]], na.rm = TRUE),
		max(data[[feature]], na.rm = TRUE),
		length.out = n_grid
	)
	
	purrr::map_dfr(grid_vals, function(val) {
		pred_data            <- data
		pred_data[[feature]] <- val
		tibble::tibble(
			standage = val,
			yhat     = mean(
				predict(model, data = pred_data)$predictions,
				na.rm = TRUE
			)
		)
	})
}

# Recode group labels (from quantile label strings) to clean "high"/"low".
# Handles both 25% and 10% sensitivity labels.
recode_group <- function(g) {
	dplyr::case_when(
		stringr::str_detect(g, regex("upper|high|warm|water.retentive",
																 ignore_case = TRUE)) ~ "high",
		stringr::str_detect(g, regex("lower|low|cold|sandy",
																 ignore_case = TRUE))  ~ "low",
		TRUE ~ NA_character_
	)
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Predictability helpers (VEcv)
# ══════════════════════════════════════════════════════════════════════════════

# VEcv: Variance Explained by Cross-validation.
# Equivalent to a cross-validated R²; scale-free and comparable across traits.
# Formula: 1 - MSE / Var(obs)
# Returns NA (not -Inf) when variance of obs is zero.
VEcv <- function(obs, pred) {
	keep <- !is.na(obs) & !is.na(pred)
	obs  <- obs[keep]
	pred <- pred[keep]
	v <- var(obs)
	if (is.na(v) || v == 0) return(NA_real_)
	1 - mean((obs - pred)^2) / v
}

# E1: Legates-McCabe efficiency statistic.
# Less sensitive to outliers than VEcv because it uses absolute rather than
# squared errors. Useful as a robustness check alongside VEcv.
# Formula: 1 - sum(|obs - pred|) / sum(|obs - mean(obs)|)
E1 <- function(obs, pred) {
	keep  <- !is.na(obs) & !is.na(pred)
	obs   <- obs[keep]
	pred  <- pred[keep]
	denom <- sum(abs(obs - mean(obs)))
	if (is.na(denom) || denom == 0) return(NA_real_)
	1 - sum(abs(obs - pred)) / denom
}

# Compute out-of-fold VEcv and E1 skill scores per stand age bin and
# environmental stratum using repeated k-fold cross-validation.
#
# For each repeat:
#   1. Create v folds from data
#   2. For each fold, train RF on the other folds, predict on held-out fold
#   3. Collect all out-of-fold predictions
#   4. Stratify by environmental variable quantiles
#   5. Compute VEcv and E1 per stratum × stand age bin
#
# Arguments:
#   trait           : character, name of response variable
#   data            : data frame containing trait + covariates
#   covariates      : character vector of predictor names (explicit — do not
#                     use setdiff(names(data), traits) which includes metadata)
#   env_vars        : character vector of environmental variables to stratify by
#   hyper_grid      : data frame with columns trait, num_trees, mtry, min_node_size
#   standage_breaks : numeric vector of bin boundaries
#   probs           : quantile probabilities for high/low stratification
#   v               : number of CV folds
#   repeats         : number of CV repeats (more = smoother uncertainty bands)
#
# Returns: long tibble with columns:
#   trait, variable, env_group, standage_bin, VEcv, E1, n, repeat_id
oof_skill_by_bins <- function(trait,
															data,
															covariates,
															env_vars,
															hyper_grid,
															standage_breaks = seq(0, 150, 10),
															probs           = c(0.25, 0.75),
															v               = 10L,
															repeats         = 30L) {
	
	stopifnot(
		"trait not in data"      = trait %in% names(data),
		"covariates not in data" = all(covariates %in% names(data)),
		"env_vars not in data"   = all(env_vars %in% names(data))
	)
	
	# Add stand age bins once (same for all repeats)
	data <- data %>%
		dplyr::mutate(
			standage_bin = cut(standage,
												 breaks         = standage_breaks,
												 include.lowest = TRUE,
												 right          = FALSE)
		)
	
	purrr::map_dfr(seq_len(repeats), function(rp) {
		
		set.seed(rp)
		
		# Create v-fold CV split indices
		fold_ids <- sample(rep(seq_len(v), length.out = nrow(data)))
		
		# Collect out-of-fold predictions across all folds
		oof_list <- vector("list", v)
		
		for (fold in seq_len(v)) {
			
			train_df <- data[fold_ids != fold, ]
			test_df  <- data[fold_ids == fold, ]
			
			# Drop rows with missing trait or covariate values in training set
			train_df <- train_df %>%
				dplyr::select(dplyr::all_of(c(trait, covariates))) %>%
				tidyr::drop_na()
			
			# Fit RF using pre-tuned hyperparameters
			rf_fit <- fit_rf_model(
				trait            = trait,
				df_train         = train_df,
				covariates       = covariates,
				hyper_parameters = hyper_grid,
				num_threads      = 1L          # no nested parallelism
			)$trait_mod
			
			# Predict on held-out fold
			test_df$oof_pred <- predict(rf_fit, data = test_df)$predictions
			oof_list[[fold]] <- test_df
		}
		
		oof <- dplyr::bind_rows(oof_list)
		
		# For each environmental variable, stratify and compute skill per bin
		purrr::map_dfr(env_vars, function(env) {
			
			qs <- quantile(oof[[env]], probs = probs, na.rm = TRUE)
			
			oof %>%
				dplyr::mutate(
					env_group = dplyr::case_when(
						.data[[env]] <= qs[1] ~ "low",
						.data[[env]] >= qs[2] ~ "high",
						TRUE ~ NA_character_
					),
					variable = env
				) %>%
				dplyr::filter(!is.na(env_group), !is.na(standage_bin)) %>%
				dplyr::group_by(
					trait     = .env$trait,
					variable,
					env_group,
					standage_bin
				) %>%
				dplyr::summarise(
					VEcv      = VEcv(.data[[trait]], oof_pred),
					E1        = E1(.data[[trait]], oof_pred),
					n         = dplyr::n(),
					.groups   = "drop"
				) %>%
				dplyr::mutate(repeat_id = rp)
		})
	})
}
