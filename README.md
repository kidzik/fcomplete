# Longitudinal data analysis using matrix completion

Suppose we observe N subjects, each subject at multiple timepoints and we want to estimate a trajectory of progression of measurements in individual subjects. For example, suppose you observe BMI of N children at different ages, as presented below

<p align="center">
   <img src="https://s3-eu-west-1.amazonaws.com/kidzinski/kidzinski/fcomplete/grouped.png" width=450 />
</div>

Here, the connected dots come from individual subjects and the black thick line corresponds to the population mean.

In this package we follow the methodology from [Kidziński, Hastie (2018)](https://arxiv.org/abs/1809.08771) to fit trajectories using matrix completion. To this end, we discretize the time grid some continous basis and find a low-rank decomposition of the dense matrix.

![Matrix completion and sparse longitudinal completion](https://s3-eu-west-1.amazonaws.com/kidzinski/kidzinski/fcomplete/intro-1.png)

In the classical matrix completion, we look for matrices `W` and `A` that fit the observed points in `Y` (green points in the image above). In our method, in order to impose smoothness, we additionaly assume the basis `B` and again we look for the reprezentation minimizing the errror. 

The interface of the package is based on the mixed-effect models in `R`. In particular, if we are given temporal observations `Y` in the long format with columns `id`, `age` and `bmi`, while additional covariates `X`, constant over time are given as a data frame with columns `id` and, say, `gender`, we can fit the model by writing

```R
model = fregression(bmi ~ age + gender | id, data = Y, covariates = X)
print(model)
```

For more information, please refer to the manual and to [vignettes](https://github.com/kidzik/fcomplete/tree/master/vignettes).

## Installation

```R
library("devtools")
install_git("https://github.com/kidzik/fcomplete/")
```

See vignettes (https://github.com/kidzik/fcomplete/tree/master/vignettes) for example use.
1. [Introduction](https://github.com/kidzik/fcomplete/tree/master/vignettes/fcomplete.ipynb)
2. [Regularization](https://github.com/kidzik/fcomplete/tree/master/vignettes/Regularization.ipynb)
3. [Regression](https://github.com/kidzik/fcomplete/tree/master/vignettes/Regression.ipynb)

## For developers

For now it's designed for editing in R Studio.

### Installation
1. Open RStudio
2. Open `fcomplete.Rproj` project
3. Install packages: `devtools`, `ggplot2`, `mclust`
4. Open the `examples/example.R` script
5. Run the script. The script will:
    * Install our `fcomplete` package
    * Load simulated data
    * Run the method
    * Plot estimates, components and projections,
    * Run basic clustering on projections

### Compiling manual

You can generate the manual by running
```bash
R CMD Rd2pdf fcomplete
```
from the parent directory to `fcomplete`.

### Running Tests

Examples are in tests for now. The main two examples (from the paper) are:
* `simulation.full.R` -- the simulation study from the paper
* `data.full.R` -- the data study from the paper
