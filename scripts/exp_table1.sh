#!/bin/bash

POLY_MODULI=(8192 16384 32768)
KEYWORD_COUNTS=(128 256 512)
HAMMING_WEIGHT=4
ITEM_LENGTH=2
OUTPUT_DIR="table1_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table1_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "N,n,m,Database Size (KB),Preparation Time (s),Total Server Time (s),Upload Cost (KB),Download Cost (KB)" > "$RESULT_FILE"

for N in "${POLY_MODULI[@]}"; do
    for n in "${KEYWORD_COUNTS[@]}"; do
        m=$N  # Set m equal to N as per experiment design
        
        LOG_FILE="$LOG_DIR/N${N}_n${n}_m${m}.log"
        
        echo "Running experiment: N=$N, n=$n, m=$m"
        echo -e "1\n0" | ./bin/bkpir -N $N -n $n -m $m -il $ITEM_LENGTH -k $HAMMING_WEIGHT > "$LOG_FILE"
        
        db_size=$(grep "Database Size:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB/*1024/; s/KB//' | bc)
        prep_time=$(grep "Preparation time:" "$LOG_FILE" | awk '{print $3}')
        server_time=$(grep "Total server time:" "$LOG_FILE" | awk '{print $4}')
        upload_cost=$(grep "Upload cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB/*1024/; s/KB//' | bc)
        download_cost=$(grep "Download cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB/*1024/; s/KB//' | bc)
        
        echo "$N,$n,$m,$db_size,$prep_time,$server_time,$upload_cost,$download_cost" >> "$RESULT_FILE"
        
        echo "Completed: N=$N, n=$n, m=$m"
    done
done

echo "All experiments completed. Results saved to $RESULT_FILE"
