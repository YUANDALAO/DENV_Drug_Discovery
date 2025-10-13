#!/bin/bash
# 完整的Git工作流：清理 + 提交
# 执行每一步都有确认和检查

set -e

echo "========================================================================"
echo "🔧 REINVENT4 项目 Git 工作流"
echo "========================================================================"

# ============================================================
# Step 1: Git配置检查
# ============================================================
echo ""
echo "📋 Step 1: 检查Git配置"
echo "----------------------------------------"

# 检查是否已初始化
if [ ! -d .git ]; then
    echo "⚠️  Git未初始化，正在初始化..."
    git init
    echo "✓ Git仓库已初始化"
else
    echo "✓ Git仓库已存在"
fi

# 检查用户配置
GIT_USER=$(git config user.name 2>/dev/null || echo "")
GIT_EMAIL=$(git config user.email 2>/dev/null || echo "")

if [ -z "$GIT_USER" ] || [ -z "$GIT_EMAIL" ]; then
    echo "⚠️  需要配置Git用户信息"
    echo "请运行:"
    echo "  git config user.name \"Your Name\""
    echo "  git config user.email \"your.email@example.com\""
    exit 1
else
    echo "✓ Git用户: $GIT_USER <$GIT_EMAIL>"
fi

# ============================================================
# Step 2: 创建.gitignore
# ============================================================
echo ""
echo "📝 Step 2: 配置.gitignore"
echo "----------------------------------------"

if [ -f .gitignore ]; then
    echo "⚠️  .gitignore已存在，备份为.gitignore.backup"
    cp .gitignore .gitignore.backup
fi

cat > .gitignore << 'GITIGNORE_EOF'
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
*.egg-info/
dist/
build/

# Jupyter
.ipynb_checkpoints

# IDE
.vscode/
.idea/
*.swp

# REINVENT4大文件
experiments/runs/*/results_*.csv
experiments/runs/*/*.csv
experiments/runs/*/*.chkpt
experiments/runs/*/*.ckpt
*.chkpt
*.ckpt
priors/*.model
priors/*.prior
!priors/.gitkeep
experiments/runs/*/tensorboard/
tb_*/
tensorboard/
*.tmp
*.log
*.out

# 保留重要文件
!configs/**/*.toml
!experiments/runs/*/config.toml
!scripts/**/*.py
!*.py
!README.md
!**/README.md
!experiments/runs/*/candidates_gold.csv
!experiments/runs/*/promising_candidates.csv
!experiments/runs/**/*.png

# Windows
Thumbs.db
Desktop.ini
$RECYCLE.BIN/
GITIGNORE_EOF

echo "✓ .gitignore已创建/更新"

# ============================================================
# Step 3: 执行项目清理
# ============================================================
echo ""
echo "🧹 Step 3: 清理项目结构"
echo "----------------------------------------"
echo "即将移动47个脚本到分类目录"
read -p "确认执行清理? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ 已取消"
    exit 1
fi

# 创建目录
mkdir -p scripts/{analysis,preprocessing,validation,utils,archive_old_scripts}
mkdir -p configs/{templates,archive}

# 核心分析脚本 (3个)
echo "  📊 分析脚本..."
[ -f unified_candidate_analysis.py ] && mv unified_candidate_analysis.py scripts/analysis/
[ -f visualize_top_structures_fixed.py ] && mv visualize_top_structures_fixed.py scripts/analysis/visualize_top_structures.py
[ -f visualize_full_analysis_fixed.py ] && mv visualize_full_analysis_fixed.py scripts/analysis/visualize_full_analysis.py

# 数据预处理
echo "  🔧 预处理脚本..."
mv clean_smiles_*.py scripts/preprocessing/ 2>/dev/null || true
mv create_*.py scripts/preprocessing/ 2>/dev/null || true
mv expand_dataset.py scripts/preprocessing/ 2>/dev/null || true
mv get_chembl_data.py scripts/preprocessing/ 2>/dev/null || true
mv prepare_*.py scripts/preprocessing/ 2>/dev/null || true
mv generate_scaffold_*.py scripts/preprocessing/ 2>/dev/null || true
mv filter_by_scaffold.py scripts/preprocessing/ 2>/dev/null || true
mv split_tl_dataset.py scripts/preprocessing/ 2>/dev/null || true

# 验证脚本
echo "  ✅ 验证脚本..."
mv validate_*.py scripts/validation/ 2>/dev/null || true
mv test_compatibility.py scripts/validation/ 2>/dev/null || true
mv correct_test.py scripts/validation/ 2>/dev/null || true
mv simple_test.py scripts/validation/ 2>/dev/null || true
mv check_scaffolds.py scripts/validation/ 2>/dev/null || true
mv analyze_scaffold_molecules.py scripts/validation/ 2>/dev/null || true

# 工具脚本
echo "  🔨 工具脚本..."
[ -f absolutely_final_smarts.py ] && mv absolutely_final_smarts.py scripts/utils/fix_smarts_final.py
mv diagnose_checkpoint.py scripts/utils/ 2>/dev/null || true
mv recover_checkpoint.py scripts/utils/ 2>/dev/null || true
mv rescore_molecules.py scripts/utils/ 2>/dev/null || true

# 归档旧脚本
echo "  📦 归档旧脚本..."
mv analyze_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv compare_runs.py scripts/archive_old_scripts/ 2>/dev/null || true
mv visualize_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv debug_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv final_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv fix_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv recursive_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv test_smarts_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv understand_*.py scripts/archive_old_scripts/ 2>/dev/null || true
mv install.py scripts/archive_old_scripts/ 2>/dev/null || true

echo "✓ 项目清理完成"

# ============================================================
# Step 4: 查看变更
# ============================================================
echo ""
echo "📊 Step 4: 查看Git变更"
echo "----------------------------------------"

git status

echo ""
echo "📈 统计:"
echo "  analysis:       $(ls -1 scripts/analysis/*.py 2>/dev/null | wc -l) 个"
echo "  preprocessing:  $(ls -1 scripts/preprocessing/*.py 2>/dev/null | wc -l) 个"
echo "  validation:     $(ls -1 scripts/validation/*.py 2>/dev/null | wc -l) 个"
echo "  utils:          $(ls -1 scripts/utils/*.py 2>/dev/null | wc -l) 个"
echo "  archive:        $(ls -1 scripts/archive_old_scripts/*.py 2>/dev/null | wc -l) 个"

# ============================================================
# Step 5: 提交到Git
# ============================================================
echo ""
echo "💾 Step 5: 提交变更"
echo "----------------------------------------"
read -p "确认提交? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ 已取消"
    exit 1
fi

# 添加所有变更
git add .

# 提交
git commit -m "refactor: Reorganize project structure

- Reorganize 47 Python scripts into categorized directories
- Add 3 core analysis scripts to scripts/analysis/
- Move preprocessing scripts to scripts/preprocessing/
- Move validation scripts to scripts/validation/
- Move utility scripts to scripts/utils/
- Archive old/deprecated scripts to scripts/archive_old_scripts/
- Configure .gitignore for large files (checkpoints, models, results)
- Improve project maintainability and clarity
"

echo "✓ 变更已提交"

# ============================================================
# Step 6: 推送到远程（可选）
# ============================================================
echo ""
echo "🚀 Step 6: 推送到远程仓库"
echo "----------------------------------------"

# 检查远程仓库
REMOTE=$(git remote -v 2>/dev/null | grep origin | head -1 || echo "")

if [ -z "$REMOTE" ]; then
    echo "⚠️  未配置远程仓库"
    echo ""
    echo "配置方法:"
    echo "  git remote add origin https://github.com/YOUR_USERNAME/REPO_NAME.git"
    echo "  git push -u origin main"
else
    echo "✓ 远程仓库: $REMOTE"
    echo ""
    read -p "是否推送到远程? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git push origin main || git push -u origin main
        echo "✓ 已推送到远程仓库"
    else
        echo "跳过推送"
    fi
fi

# ============================================================
# 完成
# ============================================================
echo ""
echo "========================================================================"
echo "✅ Git工作流完成！"
echo "========================================================================"
echo ""
echo "📂 项目结构已重组"
echo "💾 变更已提交到Git"
echo "📝 .gitignore已配置（大文件不会被提交）"
echo ""
echo "🎯 下一步:"
echo "  1. 配置远程仓库（如果还没有）"
echo "  2. 开始Run13实验配置"
echo ""
