## Defining all the necessary variables
status <- c("control", "adenoma", "crc_benign", "crc_malignant")
n_samples <- c(317, 374, 281, 252)
number_f <- c(167, 191, 139, 122)
age_mean <- c(64.52, 64.12, 65.66, 65.21)
age_sd <- c(10.12, 9.87, 10.54, 9.89)
bmi_mean <- c(25.41, 25.33, 26.01, 25.79)
bmi_sd <- c(3.78, 4.01, 3.53, 3.66)
number_smoker <- c(67, 72, 67, 66)

## Setting the seed for consistent data generation
set.seed(1597)

## Generating the data and writing it to a list
cancer <- list()

for (i in seq_along(status)) {
  prob_f <- c(number_f[i], n_samples[i] - number_f[i]) / n_samples[i]
  prob_smoke <- c(
    number_smoker[i],
    n_samples[i] - number_smoker[i]
  ) / n_samples[i]

  cancer[[i]] <- data.frame(
    status = rep(status[i], n_samples[i]),
    sex = factor(sample(0:1, n_samples[i], replace = TRUE, prob = prob_f),
      levels = c(0, 1),
      labels = c("F", "M")
    ),
    age = rnorm(n_samples[i], age_mean[i], age_sd[i]),
    bmi = rnorm(n_samples[i], bmi_mean[i], bmi_sd[i]),
    smoker = factor(
      sample(0:1, n_samples[i],
        replace = TRUE,
        prob = prob_smoke
      ),
      levels = c(0, 1),
      labels = c("yes", "no")
    )
  )
}

cancer <- Reduce(rbind, cancer)

usethis::use_data(cancer, overwrite = TRUE)
