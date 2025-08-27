#!/bin/bash

POLY_MODULUS=8192
KEYWORD_COUNT=1
ITEM_PAIRS=(
    "16384 4"
    "65536 16"
    "262144 64"
    "1048576 256"
)
HAMMING_WEIGHT=4
MODES=("MMKPIR" "AND" "OR" "NOT")
OUTPUT_DIR="table7_experiment_results"
LOG_DIR="$OUTPUT_DIR/logs"
RESULT_FILE="$OUTPUT_DIR/table7_results.csv"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "Method,m,Item Size (B),e_fl,Preparation Time (s),Selection Time (s),Logical Operation Time (s),Inner Product Time (s),Total Server Time (s),Upload Cost (KB),Download Cost (KB)" > "$RESULT_FILE"

for pair in "${ITEM_PAIRS[@]}"; do
    m=$(echo $pair | awk '{print $1}')
    il=$(echo $pair | awk '{print $2}')
    
    for mode in "${MODES[@]}"; do
        case $mode in
            "MMKPIR") menu=1 ;;
            "AND") menu=3 ;;
            "OR") menu=4 ;;
            "NOT") menu=5 ;;
        esac
        
        LOG_FILE="$LOG_DIR/m${m}_il${il}_${mode}.log"
        
        echo "Running $mode: m=$m, il=$il"
        echo -e "$menu\n0" | ./bin/bkpir -N $POLY_MODULUS -n $KEYWORD_COUNT -m $m -il $il -k $HAMMING_WEIGHT > "$LOG_FILE"
        
        prep_time=$(grep "Preparation time:" "$LOG_FILE" | awk '{print $3}')
        upload_cost=$(grep "Upload cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
        download_cost=$(grep "Download cost:" "$LOG_FILE" | awk '{print $3}' | sed 's/KB//')
        
        if [ "$mode" == "MMKPIR" ]; then
            e_fl=$(grep "e = " "$LOG_FILE" | awk '{print $3}')
        else
            line=$(grep "l=.*f =" "$LOG_FILE")
            l_value=$(echo "$line" | awk -F 'l=' '{print $2}' | awk -F ',' '{print $1}' | tr -d ' ')
            f_value=$(echo "$line" | awk -F 'f =' '{print $2}' | awk '{print $1}')
            
            if [ -z "$l_value" ]; then l_value=1; fi
            if [ -z "$f_value" ]; then f_value=1; fi
            
            e_fl=$((l_value * f_value))
        fi
        
        if [ "$mode" == "MMKPIR" ]; then
            selection_time=$(grep "Selection time:" "$LOG_FILE" | awk '{print $3}')
            logical_op_time="0.0000"  # MMKPIR doesn't have logical operations
            inner_product_time=$(grep "Inner product time:" "$LOG_FILE" | awk '{print $4}')
            total_time=$(grep "Total server time:" "$LOG_FILE" | awk '{print $4}')
        elif [ "$mode" == "NOT" ]; then
            selection_time=$(grep "Selection time:" "$LOG_FILE" | awk '{print $3}')
            logical_op_time=$(grep "NOT operation time:" "$LOG_FILE" | awk '{print $4}')
            inner_product_time=$(grep "Inner product time:" "$LOG_FILE" | awk '{print $4}')
            total_time=$(grep "BKPIR(L)_NOT total server time:" "$LOG_FILE" | awk '{print $5}')
        else
            selection_time=$(grep "Total selection time:" "$LOG_FILE" | awk '{print $4}')
            if [ "$mode" == "AND" ]; then
                logical_op_time=$(grep "AND operation time:" "$LOG_FILE" | awk '{print $4}')
            else
                logical_op_time=$(grep "OR operation time:" "$LOG_FILE" | awk '{print $4}')
            fi
            inner_product_time=$(grep "Inner product time:" "$LOG_FILE" | awk '{print $4}')
            if [ "$mode" == "AND" ]; then
                total_time=$(grep "BKPIR(L)_AND total server time:" "$LOG_FILE" | awk '{print $5}')
            else
                total_time=$(grep "BKPIR(L)_OR total server time:" "$LOG_FILE" | awk '{print $5}')
            fi
        fi
        
        echo "$mode,$m,$il,$e_fl,$prep_time,$selection_time,$logical_op_time,$inner_product_time,$total_time,$upload_cost,$download_cost" >> "$RESULT_FILE"
        
        echo "Completed: m=$m, il=$il, mode=$mode"
    done
done

echo "All Table 7 experiments completed. Raw results saved to $RESULT_FILE"
