# `fcomplete` - Functional matrix completion

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
