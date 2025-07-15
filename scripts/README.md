# BKPIR Experiment Scripts Usage Guide

This directory contains scripts for running experiments from Table 1 to Table 7. Below are instructions for running the experiments and viewing results.

## Running Experiments

1. Navigate to the build directory:
   ```bash
   cd /path/to/bkpir/build
   ```

2. Copy the scripts to the build directory:
   ```bash
   cp ../scripts/exp_table*.sh .
   ```

3. Make the scripts executable:
   ```bash
   chmod +x exp_table*.sh
   ```

4. Run each experiment script:
   ```bash
   # Run Table 1 experiments
   ./exp_table1.sh
   
   # Run Table 2 experiments
   ./exp_table2.sh
   
   # ... continue for other tables
   ./exp_table7.sh
   ```

## Viewing Results

After each experiment completes, results are saved in CSV format in their respective directories:

- Table 1: `table1_experiment_results/table1_results.csv`
- Table 2: `table2_experiment_results/table2_results.csv`
- Table 3: `table3_experiment_results/table3_results.csv`
- Table 4: `table4_experiment_results/table4_results.csv`
- Table 5: `table5_experiment_results/table5_results.csv`
- Table 6: `table6_experiment_results/table6_results.csv`
- Table 7: `table7_experiment_results/table7_results.csv`

You can view the results directly using any text editor or spreadsheet application:

```bash
# Example: View Table 1 results
cat table1_experiment_results/table1_results.csv
```

## Log Files

Detailed execution logs are available in each table's `logs/` directory:

- Table 1 logs: `table1_experiment_results/logs/`
- Table 2 logs: `table2_experiment_results/logs/`
- Table 3 logs: `table3_experiment_results/logs/`
- Table 4 logs: `table4_experiment_results/logs/`
- Table 5 logs: `table5_experiment_results/logs/`
- Table 6 logs: `table6_experiment_results/logs/`
- Table 7 logs: `table7_experiment_results/logs/`

These logs contain detailed execution information.

## Important Notes

1. **Resource Requirements**:
   - Larger experiments require significant RAM (64GB+ recommended)

2. **Execution Time**:
   - Smaller experiments take seconds to minutes
   - Larger experiments (e.g., Table 5 with larger N and k) may take hours
   - Monitor progress in the terminal output
