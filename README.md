# `fcomplete` - Functional matrix completion

For now it's designed for R Studio.

## Installation

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

## Manual

You can generate the manual by running
```bash
R CMD Rd2pdf fcomplete
```
from the parent directory to `fcomplete`.

## Tests

Examples are in tests for now. The main two examples (from the paper) are:
* `simulation.full.R` -- the simulation study from the paper
* `data.full.R` -- the data study from the paper
