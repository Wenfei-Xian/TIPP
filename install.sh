#!/bin/bash

# Exit on any error
set -e

# Initialize and update the main submodule
git submodule update --init

# Define variables for directories
KMC_DIR="src/kmc3"

# Enter the KMC directory, initialize and update its submodule, and make the project
cd "$KMC_DIR"
git submodule update --init
make

# Go back to the src directory
cd ../..

# Define variables for directories
SEQTK_DIR="src/seqtk"
cd "$SEQTK_DIR"
make

# Go back to the src directory
cd ..

# Compile the readskmercount program
g++ -o readskmercount -I./kmc3 readskmercount.opt.cpp -L./kmc3/bin -lkmc_core -pthread

chmod +x src/*.pl

echo "Installation completed successfully."
