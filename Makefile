# Makefile

# Default solver
SOLVER = Original
# Original - Original Laguerre
# ModifiedLaguerre13 - 2013 modification
# ModifiedLaguerre18 - 2018 modification

# Number type
NUMBER = double

# Exponent and mantissa used for root generator with degree greater than 4
EXPONENT = 2
MANTISSA = 4

# Degree of polynomials
DEGREE = 10

#Number of tests
N_TESTS = 1000

# Generator params. If all N_*_ROOTS are zero, then generator will create only simple roots
N_PAIRS_OF_COMPLEX_ROOTS = 0
N_CLUSTERED_ROOTS = DEGREE
N_MULTIPLE_ROOTS = 0
# Used only for clustered roots
MAX_DISTANCE_BETWEEN_CLUSTERED = 1e-5
ROOT_SWEEP_LOW = -1.0
ROOT_SWEEP_HIGH = 1.0

# For how many iterations Laguerre should trying to solve polynomial
TRIES = 80

# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall -std=c++20

# Source files
SRC = main.cpp

# Output binary
BIN = main

# Phony targets
.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -DSOLVER=$(SOLVER) -DNUMBER=$(NUMBER) -DEXPONENT=$(EXPONENT) \
	-DMANTISSA=$(MANTISSA) -DDEGREE=$(DEGREE) -DN_TESTS=$(N_TESTS) \
	-DN_PAIRS_OF_COMPLEX_ROOTS=$(N_PAIRS_OF_COMPLEX_ROOTS) -DN_CLUSTERED_ROOTS=$(N_CLUSTERED_ROOTS) \
	-DN_MULTIPLE_ROOTS=$(N_MULTIPLE_ROOTS) -DMAX_DISTANCE_BETWEEN_CLUSTERED=$(MAX_DISTANCE_BETWEEN_CLUSTERED) \
	-DROOT_SWEEP_LOW=$(ROOT_SWEEP_LOW) -DROOT_SWEEP_HIGH=$(ROOT_SWEEP_HIGH) -DTRIES=$(TRIES) -o $(BIN) $(SRC)

clean:
	rm -f $(BIN)
