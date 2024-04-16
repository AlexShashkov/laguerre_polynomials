# Laguerre Methods

This project uses various Laguerre methods to solve polynomials.

Научный руководитель: Парфенов Денис Васильевич (promasterden@yandex.ru)

КМБО-03-20: Александр Шашков (cracecrunch@gmail.com), Вячеслав Иванов (noviacheslav1@gmail.com), Глеб Дорошенко (dgleb02@mail.ru), Дмитрий Балашов (dimabalash0v@yandex.ru), Роман Охотников (roma.ohotnikov@gmail.com), Шахин Гаджиев (shakhingadzhiev@yandex.ru)

[Documentation](https://alexshashkov.github.io/laguerre_polynomials/html/md_README.html)

## MakeFile Variables

- `SOLVER`: The solver to use. Default is `2`. Options are:

  - `0`: Original Laguerre implementation
  - `1`: [On a family of Laguerre methods to find multiple roots of nonlinear equations](https://www.sciencedirect.com/science/article/abs/pii/S0096300313004918)
  - `2`: [An effective implementation of a modified Laguerre method for the roots of a polynomial](https://www.researchgate.net/publication/329607299_An_effective_implementation_of_a_modified_Laguerre_method_for_the_roots_of_a_polynomial)
- `NUMBER`: The number type to use. Default is `double`.
- `EXPONENT` and `MANTISSA`: Used for root generator. Default values are `4` and `32` respectively.
- `DEGREE`: The degree of the polynomials. Default is `100`.
- `N_TESTS`: The number of tests to run. Default is `10000`.
- `N_PAIRS_OF_COMPLEX_ROOTS`, `N_CLUSTERED_ROOTS`, `N_MULTIPLE_ROOTS`: Parameters for the generator. If all are zero, then the generator will create only simple roots. Default values are `0`, `0`, and `0` respectively.
- `MAX_DISTANCE_BETWEEN_CLUSTERED`: Used only for clustered roots. Default is `1e-5`.
- `ROOT_SWEEP_LOW` and `ROOT_SWEEP_HIGH`: The range for the root sweep. Default values are `-1.0` and `1.0` respectively.
- `TRIES`: The number of iterations Laguerre should try to solve the polynomial. Default is `80`.

## Usage

To compile the program, run `make`.

To clean up the compiled files, run `make clean`.

For example, to run modified version of Laguerre with float type for polynomials degree of 5: `make SOLVER=2 NUMBER=float DEGREE=5`

We dont recommend to use float type for polynomial generation with degree higher than 5 - resulting coefficients will loose huge amount of precision during convertation from bignum to float. Double fails to work on degree greater than or equal to 200.