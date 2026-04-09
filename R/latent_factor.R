#' Estimate latent factors using VAE
#'
#' @param E residual matrix (gene x sample)
#' @param k number of latent factors
#' @param epochs training epochs
#' @param batch_size batch size
#' @param beta KL divergence weight
#' @param seed random seed
#'
#' @return latent factor matrix (sample x k)
#' @export
estimate_latent_factors <- function(E, k = 2, epochs = 150, batch_size = 32, beta = 0.1, seed = 123){

  set.seed(seed)

  X <- t(E)
  keep <- apply(X, 2, var) > 0
  X <- X[, keep]
  X <- scale(X)

  input_dim <- ncol(X)
  latent_dim <- as.integer(k)
  intermediate_dim <- 128L

  encoder_inputs <- keras3::layer_input(shape = input_dim)

  x <- encoder_inputs %>%
    keras3::layer_dense(units = intermediate_dim, activation = "relu")

  z_mean <- x %>% keras3::layer_dense(units = latent_dim)
  z_log_var <- x %>%
    keras3::layer_dense(units = latent_dim) %>%
    keras3::layer_lambda(function(x) tensorflow::tf$clip_by_value(x, -10, 10))

  encoder <- keras3::keras_model(
    encoder_inputs,
    list(z_mean, z_log_var)
  )

  layer_sampler <- keras3::new_layer_class(
    classname = "Sampler",
    call = function(z_mean, z_log_var) {
      epsilon <- tensorflow::tf$random$normal(shape = tensorflow::tf$shape(z_mean))
      z_mean + tensorflow::tf$exp(0.5 * z_log_var) * epsilon
    }
  )

  latent_inputs <- keras3::layer_input(shape = latent_dim)

  decoder_outputs <- latent_inputs %>%
    keras3::layer_dense(units = intermediate_dim, activation = "relu") %>%
    keras3::layer_dense(units = input_dim)

  decoder <- keras3::keras_model(latent_inputs, decoder_outputs)

  model_vae <- keras3::new_model_class(
    classname = "VAE",

    initialize = function(encoder, decoder, ...) {
      super$initialize(...)
      self$encoder <- encoder
      self$decoder <- decoder
      self$sampler <- layer_sampler()

      self$total_loss_tracker <- keras3::metric_mean()
    },

    train_step = function(data) {
      with(tensorflow::tf$GradientTape() %as% tape, {

        encoded <- self$encoder(data)
        z_mean <- encoded[[1]]
        z_log_var <- encoded[[2]]

        z <- self$sampler(z_mean, z_log_var)
        reconstruction <- self$decoder(z)

        recon_loss <- tensorflow::tf$reduce_mean(tensorflow::tf$square(data - reconstruction))

        kl_loss <- -0.5 * tensorflow::tf$reduce_mean(
          1 + z_log_var - tensorflow::tf$square(z_mean) - tensorflow::tf$exp(z_log_var)
        )

        total_loss <- recon_loss + beta * kl_loss
      })

      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(
        purrr::transpose(list(grads, self$trainable_weights))
      )

      self$total_loss_tracker$update_state(total_loss)

      list(loss = self$total_loss_tracker$result())
    }
  )

  vae <- model_vae(encoder, decoder)

  vae %>% keras3::compile(
    optimizer = keras3::optimizer_adam(learning_rate = 1e-4)
  )

  vae %>% keras3::fit(
    X,
    epochs = epochs,
    batch_size = batch_size,
    verbose = 0
  )

  encoded <- predict(encoder, X)
  Z <- encoded[[1]]

  colnames(Z) <- paste0("LF", 1:k)

  Z
}
