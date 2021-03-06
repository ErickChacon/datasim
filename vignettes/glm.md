
# My doc

## Introduction

In this vignette, we show how to simulate datasets to fit generalized linear models.
We will start with a logistic linear model and then we will extend to common models
used in the GLM framework.

## Required packages

```{r}
library(datasim)
```

## Logistic Models

To simulate datasets for GLM, we need to identify the generator that can be
used to simulate the response variable and its parameters. In the case of a
logistic model, the generator is the `rbinom` function with parameters `size` and
`prob`. We need to define a list of formulas specifying the type of effect that are
included in the linear predictor of each parameter. For example, this list can be
defined as follows.

```{r}
f <- list(
  prob ~ I(1 * x1) + fa(sex, beta = c(0, -1), levels = c("male", "female")),
  size ~ I(1)
  )
```

In this formula, a linear effect on `x1` and a factor effect on `sex` are
included on the linear predictor of the `prob` parameter, while the `size` is equal
to 1 to obtain a Bernoulli density.

The data for our model can be simulated with the function `sim_model`, which requires
the inverse link function for each parameter through the argument `link_inv`, the
`generator` and the sample size `n`. In the case of a logistic model model,
the inverse link function for `prob` is the `logistic` function and identity
function for `size`.

```{r}
data_model <- sim_model(formula = f,
                        link_inv = list(psych::logistic, identity),
                        generator = rbinom,
                        n = 100, seed = 1)
```

The first 10 rows of the generated dataset looks as follows:
```{r, echo = FALSE, results = 'asis'}
knitr::kable(head(data_model, 10))
```

it contains an unique `id` for each individual, all the predictors included in the
`formula` (i.e. `x1` and `sex`), the parameters (`prob` and `size`), and the simulated
`response` variable.

Once the data is simulated, we can use it to compare models, effects, etc. For
example, we can fit a logistic model to our simulated data with:

```{r, message = FALSE}
logis <- glm(response ~ x1 + sex, binomial, data_model)
logis_sum <- confint(logis)
```
The confidence intervals can be compared with the real parameters.

```{r, echo = FALSE, results = 'asis'}
knitr::kable(logis_sum)
```

