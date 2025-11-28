#!/bin/bash
# REINVENT4 DGX Spark Quick Start Script

set -e

echo "=================================================================="
echo "REINVENT4 on NVIDIA DGX Spark - Quick Start"
echo "=================================================================="
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 设置变量
REINVENT_HOME="$HOME/aidd/REINVENT4"
EXP_NAME="spark_run1"
EXP_DIR="$REINVENT_HOME/experiments/runs/$EXP_NAME"

# 函数：打印带颜色的消息
print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_info() {
    echo -e "${YELLOW}ℹ️  $1${NC}"
}

# 检查REINVENT目录
echo "步骤 1/6: 检查REINVENT4安装..."
if [ ! -d "$REINVENT_HOME" ]; then
    print_error "REINVENT4目录不存在: $REINVENT_HOME"
    exit 1
fi
print_success "REINVENT4目录存在"
cd "$REINVENT_HOME"

# 检查conda环境
echo ""
echo "步骤 2/6: 检查conda环境..."
if ! command -v conda &> /dev/null; then
    print_error "conda未安装或不在PATH中"
    exit 1
fi

# 尝试激活环境
if conda env list | grep -q "reinvent4"; then
    print_success "找到reinvent4环境"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate reinvent4
    print_success "已激活reinvent4环境"
else
    print_warning "未找到reinvent4环境，尝试使用当前环境"
fi

# 检查必需文件
echo ""
echo "步骤 3/6: 检查必需文件..."
MISSING_FILES=0

check_file() {
    if [ -f "$1" ]; then
        print_success "$2"
    else
        print_error "$2 - 文件不存在: $1"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
}

check_file "priors/libinvent.prior" "Prior模型"
check_file "priors/denv_libinvent_model.model" "Agent模型"
check_file "data/pyrrolidine_dual_aryl.smi" "Scaffold文件"
check_file "models/random_forest_champion.joblib" "QSAR模型"
check_file "reinvent/reinvent_plugins/components/comp_qsar_scorer.py" "QSAR插件"

if [ $MISSING_FILES -gt 0 ]; then
    print_error "有 $MISSING_FILES 个文件缺失，请先复制文件"
    echo ""
    print_info "请参考 INSTALLATION_GUIDE.md 完成文件复制"
    exit 1
fi

# 创建实验目录
echo ""
echo "步骤 4/6: 创建实验目录..."
mkdir -p "$EXP_DIR"
print_success "实验目录已创建: $EXP_DIR"

# 检查配置文件
echo ""
echo "步骤 5/6: 检查配置文件..."
CONFIG_FILE="$EXP_DIR/config.toml"

if [ ! -f "$CONFIG_FILE" ]; then
    print_warning "配置文件不存在，需要手动复制"
    echo ""
    print_info "请将 spark_optimized_config.toml 复制到:"
    echo "  $CONFIG_FILE"
    echo ""
    read -p "已复制配置文件? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_error "取消启动"
        exit 1
    fi
fi

if [ -f "$CONFIG_FILE" ]; then
    print_success "配置文件存在"
else
    print_error "配置文件仍然不存在"
    exit 1
fi

# 询问运行方式
echo ""
echo "步骤 6/6: 选择运行方式"
echo "=================================================================="
echo "1. 前台运行 (实时查看输出)"
echo "2. 后台运行 (nohup)"
echo "3. 仅检查，不启动"
echo "=================================================================="
read -p "请选择 (1/2/3): " choice

case $choice in
    1)
        print_info "准备前台运行..."
        echo ""
        echo "开始训练 (按 Ctrl+C 终止)..."
        echo "=================================================================="
        reinvent "$CONFIG_FILE" 2>&1 | tee "$EXP_DIR/training.log"
        ;;
    2)
        print_info "准备后台运行..."
        echo ""
        nohup reinvent "$CONFIG_FILE" > "$EXP_DIR/training.log" 2>&1 &
        PID=$!
        print_success "训练已在后台启动 (PID: $PID)"
        echo ""
        print_info "查看日志: tail -f $EXP_DIR/training.log"
        print_info "查看进程: ps aux | grep $PID"
        print_info "终止训练: kill $PID"
        echo ""
        print_info "启动TensorBoard监控:"
        echo "  cd $REINVENT_HOME"
        echo "  tensorboard --logdir $EXP_DIR/tensorboard --bind_all"
        ;;
    3)
        print_success "环境检查完成，未启动训练"
        echo ""
        print_info "要启动训练，请运行:"
        echo "  cd $REINVENT_HOME"
        echo "  reinvent $CONFIG_FILE"
        ;;
    *)
        print_error "无效选择"
        exit 1
        ;;
esac

echo ""
echo "=================================================================="
print_success "Quick Start脚本完成"
echo "=================================================================="

