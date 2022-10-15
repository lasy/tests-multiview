
# FITS ###### 



split_and_fit <- function(x_list, y, train_test_p, family = gaussian(), seed = 1){
  
  cat(seed, "\n")
  
  # check x_list feature_names
  if (colnames(x_list$x) %>% is.null()) colnames(x_list$x) <- str_c("x_",1:ncol(x_list$x))
  if (colnames(x_list$z) %>% is.null()) colnames(x_list$z) <- str_c("z_",1:ncol(x_list$z))
  
  # split into training and test sets
  set.seed(seed)
  train_test <- create_train_test(x_list, y, train_test_p)
  x <- train_test$x_list_train$x
  z <- train_test$x_list_train$z
  
  # fit the models
  x_only <- cv.glmnet(x = x, y = train_test$y_train, family = family)
  z_only <- cv.glmnet(x = z, y = train_test$y_train, family = family)
  early <- cv.glmnet(x = cbind(x, z), y = train_test$y_train, family = family)
  late <- late_fusion(x, z, x_only, z_only, y = train_test$y_train, family = family)
  coop_0 <- cv.multiview(x_list = train_test$x_list_train, y = train_test$y_train, rho = 0, family = family)
  coop_0.5 <- cv.multiview(x_list = train_test$x_list_train, y = train_test$y_train, rho = 0.5, family = family)
  coop_1 <- cv.multiview(x_list = train_test$x_list_train, y = train_test$y_train, rho = 1, family = family)
  coop_2 <- cv.multiview(x_list = train_test$x_list_train, y = train_test$y_train, rho = 2, family = family)
  
  list(
    seed = seed, 
    train_test = train_test,
    models = list(
      x_only = x_only, 
      z_only = z_only,
      early = early,
      late = late,
      coop_0 = coop_0,
      coop_0.5 = coop_0.5,
      coop_1 = coop_1,
      coop_2 = coop_2
    )
  )
}


late_fusion <- function(x, z, x_only, z_only, y, family = gaussian()){
  l_x <- predict(x_only, newx = x, s = "lambda.min", type = "response") %>% as.vector()
  l_z <- predict(z_only, newx = z, s = "lambda.min", type = "response") %>% as.vector()
  late <- glm("y ~ l_x + l_z", data = tibble(l_x = l_x, l_z = l_z, y = y), family = family)
  late
}



create_train_test <- function(x_list, y, train_test_p){
  train_ind <- sample(seq_along(y), floor(train_test_p * length(y))) %>% sort()
  test_ind <- setdiff(seq_along(y), train_ind)
  
  list(
    x_list_train = list(x = x_list$x[train_ind, ], z = x_list$z[train_ind, ]),
    x_list_test = list(x = x_list$x[test_ind, ], z = x_list$z[test_ind, ]),
    y_train = y[train_ind],
    y_test = y[test_ind]
  )
}


# PREDICTIONS ######### 

# we want a tibble with columns : y, y_hat, set, method, seed


predict_y <- function(fits) {
  purrr::map_dfr(
    .x = seq_along(fits),
    .f = predict_y_for_fits_i,
    fits = fits
  ) %>% 
    format_method()
}

format_method <- function(df) {
  df %>% 
    mutate(
      method_group = 
        case_when(
          str_detect(method, "only") ~ "single view",
          method %in% c("early","late") ~ "simple multiview",
          str_detect(method, "coop") ~ "coop. multiview",
          TRUE ~ "???"
        ),
      method_print = method %>% str_replace("_"," "),
      method = method %>% factor(., levels = unique(method)),
      method_print = method_print %>% factor(., levels = unique(method_print))
    )
}

predict_y_for_fits_i <- function(fits, i){
  this_seed_fits <- fits[[i]]
  
  purrr::map_dfr(
    .x = names(this_seed_fits$models),
    .f = predict_y_for_method,
    this_seed_fits = this_seed_fits
  ) %>% 
    mutate(seed = i)
  
}

predict_y_for_method <- function(this_seed_fits, method){
  if (method %in% c("x_only", "z_only", "early")) {
    res <- predict_glmnet(this_seed_fits, method)
  } else if (grepl("coop", method)) {
    res <- predict_coop(this_seed_fits, method)
  } else if (method == "late") {
    res <- predict_late(this_seed_fits, method)
  } else stop(paste("no y prediction method for method: ", method))
  res
}


predict_glmnet <- function(this_seed_fits, method) {
  
  y_train <- this_seed_fits$train_test$y_train
  y_test <- this_seed_fits$train_test$y_test
  model <- this_seed_fits$models[[method]]
  
  if (method == "x_only") {
    x_train <- this_seed_fits$train_test$x_list_train$x
    x_test <-  this_seed_fits$train_test$x_list_test$x
  } else if (method == "z_only") {
    x_train <- this_seed_fits$train_test$x_list_train$z
    x_test <-  this_seed_fits$train_test$x_list_test$z
  } else if (method == "early") {
    x_train <- 
      cbind(
        this_seed_fits$train_test$x_list_train$x,
        this_seed_fits$train_test$x_list_train$z
      )
    x_test <- 
      cbind(
        this_seed_fits$train_test$x_list_test$x,
        this_seed_fits$train_test$x_list_test$z
      )
  }
  
  res <- 
    bind_rows(
      tibble(
        y = y_train, 
        y_hat =  predict(model, newx = x_train, s = "lambda.min", type = "response") %>% as.vector(), 
        set = "training"
      ),
      tibble(
        y = y_test, 
        y_hat = predict(model, newx = x_test, s = "lambda.min", type = "response") %>% as.vector(), 
        set = "test")
    ) %>% 
    mutate(method = method)
  
  res
}



predict_coop <- function(this_seed_fits, method) {
  
  y_train <- this_seed_fits$train_test$y_train
  y_test <- this_seed_fits$train_test$y_test
  model <- this_seed_fits$models[[method]]
  
  x_list_train <- this_seed_fits$train_test$x_list_train
  x_list_test <- this_seed_fits$train_test$x_list_test
  
  rho <- str_remove(method, "coop_") %>% as.numeric()
  
  res <- 
    bind_rows(
      tibble(
        y = y_train, 
        y_hat =  predict(model, newx = x_list_train, s = "lambda.min", type = "response") %>% as.vector(), 
        set = "training"
      ),
      tibble(
        y = y_test, 
        y_hat = predict(model, newx = x_list_test, s = "lambda.min", type = "response") %>% as.vector(), 
        set = "test")
    ) %>% 
    mutate(method = method)
  
  res
}


predict_late <- function(this_seed_fits, method) {
  
  y_x <- predict_glmnet(this_seed_fits, method = "x_only")
  y_z <- predict_glmnet(this_seed_fits, method = "z_only")
  
  model <- this_seed_fits$models[[method]]
  res <- 
    tibble(
      y = y_x$y,
      y_hat = predict(model, newdata = tibble(l_x = y_x$y_hat, l_z = y_z$y_hat), type = "response"),
      set = y_x$set,
      method = method
    )
  
  res
}


# COEFFICIENTS #######

# we want a df with columns : view name, index, value, method, seed

retrieve_non_zero_coeff <- function(fits) {
  
  feature_names <- 
    c(
      fits[[1]]$train_test$x_list_train$x %>% colnames(),
      fits[[1]]$train_test$x_list_train$z %>% colnames()
    )
  
  purrr::map_dfr(
    .x = seq_along(fits),
    .f = retrieve_non_zero_coeff_i,
    fits = fits
  ) %>% 
    mutate(
      name = name %>% factor(., levels = feature_names),
      view = str_sub(name,1,1)
      ) %>% 
    format_method()
}

retrieve_non_zero_coeff_i <- function(fits, i){
  this_seed_fits <- fits[[i]]
  
  purrr::map_dfr(
    .x = names(this_seed_fits$models),
    .f = retrieve_coeff_for_method,
    this_seed_fits = this_seed_fits
  ) %>% 
    mutate(seed = i)
  
}


retrieve_coeff_for_method <- function(this_seed_fits, method){
  res = tibble()
  if(method != "late") {
    m <- 
      coef(this_seed_fits$models[[method]], s = "lambda.min") %>% 
      as.matrix() 
    res <- 
      m %>% 
      as_tibble() %>% 
      mutate(
        name = rownames(m),
        index = row_number() - 1
      ) %>% 
      dplyr::rename(value = s1) %>% 
      filter(value != 0, index > 0)  %>% 
      select(name, index, value) %>% 
      mutate(method = method)
  }
  res
}


# VIZ #######


plot_predictions <- function(predicted) {
  
  ggplot(predicted, aes(x = y, y = y_hat, col = set)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(size = 0.5) +
    facet_grid(method_print ~ seed) +
    theme(legend.position = "bottom")
}

plot_MSE <- function(predicted){
  MSE <- 
    predicted %>% 
    group_by(method, set, seed) %>% 
    summarize(MSE = mean((y - y_hat)^2), .groups = "drop") %>% 
    left_join(predicted %>% select(starts_with("method")) %>% distinct(), by = "method")
  
  ggplot(MSE, aes(x = method_print, y = MSE, fill = method_group)) +
    geom_boxplot() +
    geom_point() +
    facet_grid(set ~ .) +
    guides(fill = "none")
} 




plot_nonzero_coeff <- function(selected_features) {
  
  ggplot(selected_features, 
         aes(x = value, y = name %>% forcats::fct_rev(), col = view)) +
    geom_vline(xintercept = 0) +
    geom_boxplot(outlier.size =  0.5, varwidth = TRUE) +
    facet_grid(view ~ method, scales = "free_y", space = "free_y") +
    guides(col = "none") +
    ylab("features")
  
}