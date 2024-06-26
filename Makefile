# Makefile

# Default solver
SOLVER = 2
# 0 - Original Laguerre
# 1 - 2013 modification
# 2 - 2018 modification

# Number type
NUMBER = double

# Exponent and mantissa used for root generator with degree greater than 4
# You can check recommended precision for ttmath on https://www.ttmath.org/online_calculator/
# lets use big precision
# FOR 64-bit 
# mantissa = 2048 bits / 64 bits per word = 32 words
# exponent = 256 bits / 64 bits per word = 4 words
# FOR 32-bit
# mantissa = 2048 bits / 32 bits per word = 64 words
# exponent = 256 bits / 32 bits per word = 8 words
EXPONENT = 4
MANTISSA = 32

# Degree of polynomials
DEGREE = 100

#Number of tests
N_TESTS = 10000

# Generator params. If all N_*_ROOTS are zero, then generator will create only simple roots
N_PAIRS_OF_COMPLEX_ROOTS = 0
N_CLUSTERED_ROOTS = 0
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
