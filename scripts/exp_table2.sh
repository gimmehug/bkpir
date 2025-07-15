#!/bin/bash

POLY_MODULI=(8192 16384 32768)
KEYWORD_COUNT=1
ITEM_COUNT=8192
HAMMING_WEIGHT=4
ITEM_LENGTH=2
OUTPUT_DIR="table2_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table2_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "N,Preparation Time (s),Total Server Time (s),Upload Cost (KB),Download Cost (KB)" > "$RESULT_FILE"

for N in "${POLY_MODULI[@]}"; do
    LOG_FILE="$LOG_DIR/N${N}.log"
    
    echo "Running experiment: N=$N, n=$KEYWORD_COUNT, m=$ITEM_COUNT"
    echo -e "1\n0" | ./bin/bkpir -N $N -n $KEYWORD_COUNT -m $ITEM_COUNT -il $ITEM_LENGTH -k $HAMMING_WEIGHT > "$LOG_FILE"
    
    prep_time=$(grep "Preparation time:" "$LOG_FILE" | awk '{print $3}')
    server_time=$(grep "Total server time:" "$LOG_FILE" | awk '{print $4}')
    upload_cost=$(grep "Upload cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB/*1024/; s/KB//' | bc)
    download_cost=$(grep "Download cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB/*1024/; s/KB//' | bc)
    
    echo "$N,$prep_time,$server_time,$upload_cost,$download_cost" >> "$RESULT_FILE"
    
    echo "Completed: N=$N"
done

echo "All experiments completed. Results saved to $RESULT_FILE"
