#!/bin/bash
# REINVENT4项目自动清理脚本
# 执行前会显示分类结果，需要你确认

set -e

echo "========================================================================"
echo "🧹 REINVENT4 项目清理脚本"
echo "========================================================================"

# 创建目录结构
echo ""
echo "📁 创建标准目录结构..."
mkdir -p scripts/analysis
mkdir -p scripts/preprocessing
mkdir -p scripts/validation
mkdir -p scripts/utils
mkdir -p scripts/archive_old_scripts
mkdir -p configs/templates
mkdir -p configs/archive

echo "✓ 目录创建完成"

# ============================================================
# 分类1: 核心分析脚本 (3个) → scripts/analysis/
# ============================================================
echo ""
echo "📊 [分析脚本] 移动到 scripts/analysis/"
echo "----------------------------------------"

# 保留的3个核心脚本
mv unified_candidate_analysis.py scripts/analysis/ 2>/dev/null && echo "  ✓ unified_candidate_analysis.py"
mv visualize_top_structures_fixed.py scripts/analysis/visualize_top_structures.py 2>/dev/null && echo "  ✓ visualize_top_structures.py (已重命名)"
mv visualize_full_analysis_fixed.py scripts/analysis/visualize_full_analysis.py 2>/dev/null && echo "  ✓ visualize_full_analysis.py (已重命名)"

# ============================================================
# 分类2: 数据预处理脚本 → scripts/preprocessing/
# ============================================================
echo ""
echo "🔧 [数据预处理] 移动到 scripts/preprocessing/"
echo "----------------------------------------"

# SMILES清理系列
mv clean_smiles_*.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ clean_smiles_*.py (5个)"

# 数据集创建/扩展
mv create_full_dataset.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ create_full_dataset.py"
mv create_ultra_clean_dataset.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ create_ultra_clean_dataset.py"
mv expand_dataset.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ expand_dataset.py"
mv get_chembl_data.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ get_chembl_data.py"

# Scaffold准备系列
mv prepare_scaffolds.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ prepare_scaffolds.py"
mv prepare_decorations.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ prepare_decorations.py"
mv generate_scaffold_library.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ generate_scaffold_library.py"
mv generate_scaffold_library_fixed.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ generate_scaffold_library_fixed.py"

# LibInvent数据准备
mv prepare_libinvent_*.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ prepare_libinvent_*.py (4个)"
mv filter_by_scaffold.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ filter_by_scaffold.py"

# Transfer Learning数据分割
mv split_tl_dataset.py scripts/preprocessing/ 2>/dev/null && echo "  ✓ split_tl_dataset.py"

# ============================================================
# 分类3: 验证脚本 → scripts/validation/
# ============================================================
echo ""
echo "✅ [验证脚本] 移动到 scripts/validation/"
echo "----------------------------------------"

# 分子验证
mv validate_*.py scripts/validation/ 2>/dev/null && echo "  ✓ validate_*.py (3个)"
mv test_compatibility.py scripts/validation/ 2>/dev/null && echo "  ✓ test_compatibility.py"
mv correct_test.py scripts/validation/ 2>/dev/null && echo "  ✓ correct_test.py"
mv simple_test.py scripts/validation/ 2>/dev/null && echo "  ✓ simple_test.py"

# Scaffold验证
mv check_scaffolds.py scripts/validation/ 2>/dev/null && echo "  ✓ check_scaffolds.py"
mv analyze_scaffold_molecules.py scripts/validation/ 2>/dev/null && echo "  ✓ analyze_scaffold_molecules.py"

# ============================================================
# 分类4: SMARTS修复工具 → scripts/utils/
# ============================================================
echo ""
echo "🔨 [工具脚本] 移动到 scripts/utils/"
echo "----------------------------------------"

# SMARTS修复系列（保留最终版本，其他归档）
mv absolutely_final_smarts.py scripts/utils/fix_smarts_final.py 2>/dev/null && echo "  ✓ fix_smarts_final.py (重命名自absolutely_final_smarts.py)"

# Checkpoint工具
mv diagnose_checkpoint.py scripts/utils/ 2>/dev/null && echo "  ✓ diagnose_checkpoint.py"
mv recover_checkpoint.py scripts/utils/ 2>/dev/null && echo "  ✓ recover_checkpoint.py"

# 重新打分工具
mv rescore_molecules.py scripts/utils/ 2>/dev/null && echo "  ✓ rescore_molecules.py"

# ============================================================
# 分类5: 旧版本/临时脚本 → scripts/archive_old_scripts/
# ============================================================
echo ""
echo "📦 [归档脚本] 移动到 scripts/archive_old_scripts/"
echo "----------------------------------------"

# 旧版分析脚本
mv analyze_results.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ analyze_results.py (旧版)"
mv analyze_run3_scaffolds.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ analyze_run3_scaffolds.py"
mv analyze_run4_results.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ analyze_run4_results.py"
mv analyze_top20.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ analyze_top20.py"
mv analyze_top_candidates.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ analyze_top_candidates.py"
mv compare_runs.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ compare_runs.py"

# 旧版可视化
mv visualize_results.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ visualize_results.py"
mv visualize_full_analysis.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ visualize_full_analysis.py (旧版)"
mv visualize_run3_scaffold_top20.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ visualize_run3_scaffold_top20.py"
mv visualize_top10.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ visualize_top10.py"
mv visualize_top10_structures.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ visualize_top10_structures.py"

# SMARTS调试/测试脚本（大量临时文件）
mv debug_smarts_step_by_step.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ debug_smarts_step_by_step.py"
mv final_correct_smarts.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ final_correct_smarts.py"
mv final_smarts_fix.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ final_smarts_fix.py"
mv fix_aromatic_n_bond.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ fix_aromatic_n_bond.py"
mv fix_smarts.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ fix_smarts.py"
mv recursive_smarts_approach.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ recursive_smarts_approach.py"
mv test_smarts_now.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ test_smarts_now.py"
mv validate_smarts_run4.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ validate_smarts_run4.py"
mv understand_real_structure.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ understand_real_structure.py"

# 安装脚本（通常不需要）
mv install.py scripts/archive_old_scripts/ 2>/dev/null && echo "  ✓ install.py"

# ============================================================
# 生成归档清单
# ============================================================
echo ""
echo "📋 生成归档清单..."
cat > scripts/archive_old_scripts/README.txt << 'EOF'
旧脚本归档目录
================

归档日期: 2025-10-13
原因: 项目清理，保留核心功能脚本

这些脚本已被更新的版本替代或不再使用，但保留备份以防需要参考。

如需使用归档脚本，请先检查是否有新版本在：
- scripts/analysis/
- scripts/preprocessing/
- scripts/validation/
- scripts/utils/
