#!/bin/bash

POLY_MODULI=(8192 16384 32768)
KEYWORD_COUNT=2
ITEM_LENGTH=2
HAMMING_WEIGHTS=(2 4 8 16)
OUTPUT_DIR="table7_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table7_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "N,Init Budget (bits),k,Sel+IP Budget (bits),Per-Op Budget (bits),Max Ops,Per-Op Latency (s)" > "$RESULT_FILE"

for N in "${POLY_MODULI[@]}"; do
    for k in "${HAMMING_WEIGHTS[@]}"; do
        LOG_FILE="$LOG_DIR/N${N}_k${k}.log"
        m=$N
        echo "Running OR: N=$N, k=$k"
        echo -e "4\n0" | ./bin/bkpir -N $N -n $KEYWORD_COUNT -m $m -il $ITEM_LENGTH -k $k > "$LOG_FILE"
        
        init_budget=$(grep "Fresh noise budget:" "$LOG_FILE" | awk '{print $4}')
        sel_ip_budget=$(grep "Basic budget (selection + inner product) cost:" "$LOG_FILE" | awk '{print $8}')
        per_op_budget=$(grep "Logic noise budget cost:" "$LOG_FILE" | awk '{print $5}')
        max_ops=$(grep "MAX logic depth supported:" "$LOG_FILE" | awk '{print $5}')
        per_op_latency=$(grep "OR operation time:" "$LOG_FILE" | awk '{print $4}')
        
        echo "$N,$init_budget,$k,$sel_ip_budget,$per_op_budget,$max_ops,$per_op_latency" >> "$RESULT_FILE"
        
        echo "Completed: N=$N, k=$k"
    done
done

echo "All Table 7 experiments completed. Raw results saved to $RESULT_FILE"
