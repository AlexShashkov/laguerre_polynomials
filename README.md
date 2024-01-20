# Laguerre Methods

This project uses various Laguerre methods to solve polynomials.

## MakeFile Variables

- `SOLVER`: The solver to use. Default is `Original`. Options are:

  - `Original`: Original Laguerre
  - `ModifiedLaguerre13`: 2013 modification
  - `ModifiedLaguerre18`: 2018 modification
- `NUMBER`: The number type to use. Default is `double`.
- `EXPONENT` and `MANTISSA`: Used for root generator with degree greater than 4. Default values are `2` and `4` respectively.
- `DEGREE`: The degree of the polynomials. Default is `10`.
- `N_TESTS`: The number of tests to run. Default is `1000`.
- `N_PAIRS_OF_COMPLEX_ROOTS`, `N_CLUSTERED_ROOTS`, `N_MULTIPLE_ROOTS`: Parameters for the generator. If all are zero, then the generator will create only simple roots. Default values are `0`, `DEGREE`, and `0` respectively.
- `MAX_DISTANCE_BETWEEN_CLUSTERED`: Used only for clustered roots. Default is `1e-5`.
- `ROOT_SWEEP_LOW` and `ROOT_SWEEP_HIGH`: The range for the root sweep. Default values are `-1.0` and `1.0` respectively.
- `TRIES`: The number of iterations Laguerre should try to solve the polynomial. Default is `80`.

## Usage

To compile the program, run `make`.

To clean up the compiled files, run `make clean`.

For example, to run modified version of Laguerre with float type of polynomials degree of 15 run `make SOLVER=ModifiedLaguerre13 NUMBER=float DEGREE=15`
