#!/bin/bash
echo "实验对比"
echo "===================="
for dir in runs/*/; do
    name=$(basename $dir)
    if [ -f "${dir}results.csv" ]; then
        echo "[$name]"
        python3 -c "
import pandas as pd
df = pd.read_csv('${dir}results.csv')
print(f'  Score: {df[\"Score\"].max():.3f}')
print(f'  pIC50: {df[\"Activity (raw)\"].mean():.2f}')
print(f'  分子数: {len(df):,}')
"
        echo ""
    fi
done
