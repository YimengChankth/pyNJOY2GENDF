#!/bin/bash

echo "=========================================="
echo "Running MatXS (main.py)"
echo "=========================================="
python3 main.py

echo ""
echo "=========================================="
echo "Running MONTE_extractmacro"
echo "=========================================="
cd output_multiple/MONTE_extractmacro_comparison
./run.sh
cd ../..

echo ""
echo "=========================================="
echo "Comparison complete"
echo "=========================================="
