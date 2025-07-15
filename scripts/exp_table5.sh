#!/bin/bash

POLY_MODULI=(8192 16384 32768)
KEYWORD_COUNT=128
ITEM_COUNT=8192
ITEM_LENGTH=2
HAMMING_WEIGHTS=(2 4 8 16)
OUTPUT_DIR="table5_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table5_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "N,k,Mode,Server Time (s)" > "$RESULT_FILE"

for N in "${POLY_MODULI[@]}"; do
    for k in "${HAMMING_WEIGHTS[@]}"; do
        if [ $N -eq 8192 ] && [ $k -gt 4 ]; then
            continue
        fi
        
        MMKPIR_LOG="$LOG_DIR/N${N}_k${k}_mmkpir.log"
        echo "Running MMKPIR: N=$N, k=$k"
        echo -e "1\n0" | ./bin/bkpir -N $N -n $KEYWORD_COUNT -m $ITEM_COUNT -il $ITEM_LENGTH -k $k > "$MMKPIR_LOG"
        mmkpir_time=$(grep "Total server time:" "$MMKPIR_LOG" | awk '{print $4}')
        echo "$N,$k,MMKPIR,$mmkpir_time" >> "$RESULT_FILE"
        
        BKPIR_LOG="$LOG_DIR/N${N}_k${k}_bkpir.log"
        echo "Running BKPIR: N=$N, k=$k"
        echo -e "2\n0" | ./bin/bkpir -N $N -n $KEYWORD_COUNT -m $ITEM_COUNT -il $ITEM_LENGTH -k $k > "$BKPIR_LOG"
        bkpir_time=$(grep "BKPIR(L) total server time:" "$BKPIR_LOG" | awk '{print $4}')
        echo "$N,$k,BKPIR,$bkpir_time" >> "$RESULT_FILE"
        
        echo "Completed: N=$N, k=$k"
    done
done

echo "All Table 4 experiments completed. Raw results saved to $RESULT_FILE"
