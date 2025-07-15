#!/bin/bash

POLY_MODULUS=8192
KEYWORD_COUNTS=(128 256 512)
ITEM_COUNTS=(8192 16384 32768 65536)
ITEM_LENGTH=2
HAMMING_WEIGHT=4
MODES=("AND" "OR" "NOT")
OUTPUT_DIR="table4_expriment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table4_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "n,m,Mode,Preparation Time (s),Selection Time (s),Logical Operation Time (s),Inner Product Time (s),Total Server Time (s),Upload Cost (KB),Download Cost (KB)" > "$RESULT_FILE"

for n in "${KEYWORD_COUNTS[@]}"; do
    for m in "${ITEM_COUNTS[@]}"; do
        for mode in "${MODES[@]}"; do
            case $mode in
                "AND") menu=3 ;;
                "OR") menu=4 ;;
                "NOT") menu=5 ;;
            esac
            
            LOG_FILE="$LOG_DIR/n${n}_m${m}_${mode}.log"
            
            echo "Running $mode: n=$n, m=$m"
            echo -e "$menu\n0" | ./bin/bkpir -N $POLY_MODULUS -n $n -m $m -il $ITEM_LENGTH -k $HAMMING_WEIGHT > "$LOG_FILE"
            
            prep_time=$(grep "Preparation time:" "$LOG_FILE" | awk '{print $3}')
            
            if [ "$mode" == "NOT" ]; then
                selection_time=$(grep "Selection time:" "$LOG_FILE" | awk '{print $3}')
            else
                selection_time=$(grep "Total selection time:" "$LOG_FILE" | awk '{print $4}')
            fi
            
            case $mode in
                "AND") logical_op_time=$(grep "AND operation time:" "$LOG_FILE" | awk '{print $4}') ;;
                "OR") logical_op_time=$(grep "OR operation time:" "$LOG_FILE" | awk '{print $4}') ;;
                "NOT") logical_op_time=$(grep "NOT operation time:" "$LOG_FILE" | awk '{print $4}') ;;
            esac
            
            inner_product_time=$(grep "Inner product time:" "$LOG_FILE" | awk '{print $4}')
            
            case $mode in
                "AND") total_time=$(grep "BKPIR(L)_AND total server time:" "$LOG_FILE" | awk '{print $5}') ;;
                "OR") total_time=$(grep "BKPIR(L)_OR total server time:" "$LOG_FILE" | awk '{print $5}') ;;
                "NOT") total_time=$(grep "BKPIR(L)_NOT total server time:" "$LOG_FILE" | awk '{print $5}') ;;
            esac
            
            upload_cost=$(grep "Upload cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
            download_cost=$(grep "Download cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
            
            echo "$n,$m,$mode,$prep_time,$selection_time,$logical_op_time,$inner_product_time,$total_time,$upload_cost,$download_cost" >> "$RESULT_FILE"
            
            echo "Completed: n=$n, m=$m, mode=$mode"
        done
    done
done

echo "All Table 5 experiments completed. Raw results saved to $RESULT_FILE"
