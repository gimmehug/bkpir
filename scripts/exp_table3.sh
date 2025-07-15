#!/bin/bash

POLY_MODULUS=8192
KEYWORD_COUNTS=(128 256 512)
ITEM_COUNT=65536
ITEM_LENGTHS=(2 4 8 16)
HAMMING_WEIGHT=4
OUTPUT_DIR="table3_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table3_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "n,m,Item Length (B),Database Size (MB),Preparation Time (s),Total Server Time (s),Upload Cost (KB),Download Cost (KB)" > "$RESULT_FILE"

for n in "${KEYWORD_COUNTS[@]}"; do
    for il in "${ITEM_LENGTHS[@]}"; do
        LOG_FILE="$LOG_DIR/n${n}_il${il}.log"
        
        echo "Running experiment: n=$n, m=$ITEM_COUNT, il=$il"
        echo -e "1\n0" | ./bin/bkpir -N $POLY_MODULUS -n $n -m $ITEM_COUNT -il $il -k $HAMMING_WEIGHT > "$LOG_FILE"
        
        db_size=$(grep "Database Size:" "$LOG_FILE" | awk '{print $3}' | sed 's/MB//')
        prep_time=$(grep "Preparation time:" "$LOG_FILE" | awk '{print $3}')
        server_time=$(grep "Total server time:" "$LOG_FILE" | awk '{print $4}')
        upload_cost=$(grep "Upload cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
        download_cost=$(grep "Download cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
        
        echo "$n,$ITEM_COUNT,$il,$db_size,$prep_time,$server_time,$upload_cost,$download_cost" >> "$RESULT_FILE"
        
        echo "Completed: n=$n, il=$il"
    done
done

echo "All Table 3 experiments completed. Results saved to $RESULT_FILE"
