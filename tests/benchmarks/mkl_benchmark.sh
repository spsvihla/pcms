#!/bin/bash
THREAD_COUNTS=(1 2 4 6 8 10 12 14 16)

for t in "${THREAD_COUNTS[@]}"; do
    echo "Running with MKL_NUM_THREADS=$t OMP_NUM_THREADS=$t"
    MKL_NUM_THREADS=$t OMP_NUM_THREADS=$t python3 mkl_benchmark.py
    echo "--------------------------------------------"
done
