#!/bin/bash
# 安全的Git提交脚本 - 自动检查大文件

set -e

echo "========================================================================"
echo "🔒 安全Git提交脚本"
echo "========================================================================"

# 完善.gitignore
echo "📝 更新.gitignore..."
cat >> .gitignore << 'EOF'

# ============================================================
# 超大文件排除（自动生成于2025-10-13）
# ============================================================

# TensorBoard（GB级）
**/events.out.tfevents.*
tb_*/
**/tensorboard*/
**/tensorboard_*/

# Checkpoints（100MB+）
**/*.chkpt
**/*.ckpt

# Prior models（77-91MB each）
priors/*.prior
priors/*.model
!priors/.gitkeep

# 大型CSV结果
experiments/runs/**/results_*.csv
experiments/runs/**/candidates_*.csv
experiments/runs/**/*_molecules*.csv
experiments/runs/**/checkpoint_*.chkpt
experiments/runs/**/diagnostic_results/
experiments/runs/**/full_evaluation/

# 主目录临时文件
denv_*.csv
backup_*.csv
*_results_*.csv
*.smi
!data/*.smi
!data/**/*.smi

# External repos
autogrow4/
contrib/
