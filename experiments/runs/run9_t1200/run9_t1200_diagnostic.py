"""
Run9_t1200 专用诊断脚本
适用域(AD)分析 - 检测生成分子是否可靠
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import os

# ============================================================================
# 配置区 - 根据您的实际情况已配置好
# ============================================================================

# 训练集文件（实际使用的清洁版本）
TRAINING_FILE = "../../../data/denv_ultra_clean.tsv"  # 如果这个文件不存在，会尝试NS3.csv

# 备选训练集（原始数据）
BACKUP_TRAINING_FILE = "../../../data/NS3.csv"

# 生成分子文件
GENERATED_FILE = "results_1.csv"

# 输出目录
OUTPUT_DIR = "diagnostic_results"

# 项目名称
PROJECT_NAME = "Run9_t1200 - DENV NS3 Inhibitor Generation"

# ============================================================================
# 数据读取函数
# ============================================================================

def load_training_data():
    """智能加载训练集数据"""
    
    print("正在加载训练集...")
    
    # 尝试加载清洁版本
    if os.path.exists(TRAINING_FILE):
        print(f"✓ 找到清洁训练集: {TRAINING_FILE}")
        try:
            # TSV文件，没有header
            with open(TRAINING_FILE, 'r') as f:
                smiles_list = [line.strip() for line in f if line.strip()]
            print(f"  加载了 {len(smiles_list)} 个SMILES")
            return smiles_list
        except Exception as e:
            print(f"  读取失败: {e}")
    
    # 尝试加载原始CSV
    print(f"尝试备选文件: {BACKUP_TRAINING_FILE}")
    try:
        df = pd.read_csv(BACKUP_TRAINING_FILE)
        
        # 尝试不同的列名
        smiles_col = None
        for col in ['Smiles', 'SMILES', 'smiles', 'canonical_smiles']:
            if col in df.columns:
                smiles_col = col
                break
        
        if smiles_col is None:
            print(f"✗ 错误: 找不到SMILES列，可用列: {df.columns.tolist()}")
            return None
        
        smiles_list = df[smiles_col].dropna().tolist()
        print(f"✓ 从 {BACKUP_TRAINING_FILE} 加载了 {len(smiles_list)} 个SMILES")
        
        # 应用超严格筛选（与您的create_ultra_clean_dataset.py一致）
        print("  应用超严格筛选...")
        clean_smiles = apply_ultra_clean_filter(smiles_list)
        print(f"  筛选后: {len(clean_smiles)} 个SMILES")
        
        return clean_smiles
        
    except Exception as e:
        print(f"✗ 读取失败: {e}")
        return None


def apply_ultra_clean_filter(smiles_list):
    """应用与训练时相同的超严格筛选"""
    
    def is_ultra_compatible(smiles):
        forbidden_tokens = ['[n-]', '[s+]', '[N-]', '[S-]', '[O+]', '[C-]', '[c-]']
        
        for token in forbidden_tokens:
            if token in smiles:
                return False
        
        if any(x in smiles for x in ['[nH+]', '[s+]', '[o+]']):
            return False
        
        if any(x in smiles for x in ['@', '/', '\\']):
            return False
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            
            if mol.GetNumAtoms() > 50:
                return False
            
            atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            allowed_atoms = {'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'H'}
            if not atoms.issubset(allowed_atoms):
                return False
            
            return True
        except:
            return False
    
    valid_smiles = []
    
    for smiles in tqdm(smiles_list, desc="筛选训练集"):
        if pd.isna(smiles):
            continue
        
        smiles_str = str(smiles).strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                if ('.' not in canonical_smiles and 
                    '|' not in canonical_smiles and
                    len(canonical_smiles) < 100 and
                    is_ultra_compatible(canonical_smiles)):
                    valid_smiles.append(canonical_smiles)
        except:
            continue
    
    # 去重
    unique_smiles = list(dict.fromkeys(valid_smiles))
    return unique_smiles


def load_generated_data():
    """加载生成的分子数据"""
    
    print(f"\n正在加载生成分子: {GENERATED_FILE}")
    
    try:
        # 明确指定逗号分隔
        df = pd.read_csv(GENERATED_FILE, sep=',')
        
        print(f"  文件形状: {df.shape}")
        print(f"  前5列: {df.columns.tolist()[:5]}")
        
        # 直接使用SMILES列
        if 'SMILES' not in df.columns:
            print(f"  ✗ 错误: 未找到SMILES列")
            print(f"  可用列: {df.columns.tolist()}")
            return None, None
        
        smiles_list = df['SMILES'].dropna().tolist()
        print(f"✓ 加载了 {len(smiles_list)} 个生成分子")
        print(f"  第1个: {smiles_list[0][:60]}...")
        
        # 读取额外信息
        extra_cols = {}
        if 'Score' in df.columns:
            extra_cols['total_score'] = df['Score'].tolist()
            print(f"  ✓ 找到Score列")
        if 'DENV_Activity (raw)' in df.columns:
            extra_cols['DENV_Activity'] = df['DENV_Activity (raw)'].tolist()
            print(f"  ✓ 找到DENV_Activity列")
        if 'NLL' in df.columns:
            extra_cols['NLL'] = df['NLL'].tolist()
        
        return smiles_list, extra_cols
        
    except Exception as e:
        print(f"✗ 读取失败: {e}")
        import traceback
        traceback.print_exc()
        return None, None


# ============================================================================
# AD分析核心函数
# ============================================================================

def compute_fingerprints(smiles_list, radius=2, n_bits=2048):
    """计算Morgan指纹"""
    fps = []
    valid_smiles = []
    
    for smi in tqdm(smiles_list, desc="计算指纹"):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            fps.append(fp)
            valid_smiles.append(smi)
    
    return fps, valid_smiles


def calculate_ad_analysis(training_fps, generated_smiles, generated_fps, extra_cols=None):
    """计算AD分析"""
    results = []
    
    print("\n计算AD距离...")
    for idx, (smi, qfp) in enumerate(tqdm(zip(generated_smiles, generated_fps), 
                                          total=len(generated_smiles))):
        # 计算到所有训练集分子的Tanimoto相似度
        similarities = [
            DataStructs.TanimotoSimilarity(qfp, tfp) 
            for tfp in training_fps
        ]
        
        max_sim = np.max(similarities)
        mean_sim = np.mean(similarities)
        
        # AD判定
        if max_sim > 0.6:
            status = "Inside (Safe)"
            risk_score = 1
        elif max_sim > 0.4:
            status = "Boundary (Caution)"
            risk_score = 3
        else:
            status = "Outside (High Risk)"
            risk_score = 5
        
        result = {
            'SMILES': smi,
            'MaxTanimoto': max_sim,
            'MeanTanimoto': mean_sim,
            'AD_Status': status,
            'Risk_Score': risk_score
        }
        
        # 添加额外信息
        if extra_cols:
            if 'total_score' in extra_cols:
                result['Total_Score'] = extra_cols['total_score'][idx]
            if 'NLL' in extra_cols:
                result['NLL'] = extra_cols['NLL'][idx]
        
        results.append(result)
    
    df = pd.DataFrame(results)
    
    # 统计报告
    print("\n" + "="*70)
    print("适用域(AD)分析报告")
    print("="*70)
    print(f"总分子数: {len(df)}")
    print(f"\nAD分布:")
    for status in df['AD_Status'].value_counts().items():
        print(f"  {status[0]}: {status[1]} ({status[1]/len(df)*100:.1f}%)")
    print(f"\n最大Tanimoto统计:")
    print(df['MaxTanimoto'].describe())
    
    if 'Total_Score' in df.columns:
        print(f"\nTotal Score统计:")
        print(df['Total_Score'].describe())
        
        # 高分但高风险的分子
        high_score_high_risk = df[(df['Total_Score'] > 0.7) & (df['Risk_Score'] >= 3)]
        print(f"\n高分但高风险分子: {len(high_score_high_risk)} 个")
    
    print("="*70)
    
    return df


def plot_ad_distribution(df, save_path):
    """可视化AD分布"""
    
    # 判断是否有score信息
    has_score = 'Total_Score' in df.columns
    
    if has_score:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
    
    # 1. Tanimoto分布直方图
    ax = axes[0]
    ax.hist(df['MaxTanimoto'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(x=0.6, color='green', linestyle='--', linewidth=2, label='Safe (>0.6)')
    ax.axvline(x=0.4, color='orange', linestyle='--', linewidth=2, label='Caution (0.4-0.6)')
    ax.set_xlabel('Max Tanimoto Similarity', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('AD Distribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    
    # 2. AD状态饼图
    ax = axes[1]
    status_counts = df['AD_Status'].value_counts()
    colors = {'Inside (Safe)': '#2ecc71', 
              'Boundary (Caution)': '#f39c12', 
              'Outside (High Risk)': '#e74c3c'}
    ax.pie(status_counts, labels=status_counts.index, autopct='%1.1f%%',
            colors=[colors.get(s, 'gray') for s in status_counts.index], 
            startangle=90, textprops={'fontsize': 9, 'fontweight': 'bold'})
    ax.set_title('AD Status', fontsize=12, fontweight='bold')
    
    # 3. 累积分布
    ax = axes[2]
    sorted_sim = np.sort(df['MaxTanimoto'])
    cumulative = np.arange(1, len(sorted_sim) + 1) / len(sorted_sim)
    ax.plot(sorted_sim, cumulative, linewidth=2, color='darkblue')
    ax.axvline(x=0.4, color='orange', linestyle='--', alpha=0.7, linewidth=2)
    ax.axvline(x=0.6, color='green', linestyle='--', alpha=0.7, linewidth=2)
    ax.set_xlabel('Max Tanimoto', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax.set_title('Cumulative Distribution', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # 4. Max vs Mean Tanimoto
    ax = axes[3]
    scatter = ax.scatter(df['MeanTanimoto'], df['MaxTanimoto'], 
                         c=df['MaxTanimoto'], cmap='RdYlGn', 
                         alpha=0.6, s=20, edgecolors='black', linewidth=0.5)
    ax.axhline(y=0.6, color='green', linestyle='--', alpha=0.5, linewidth=2)
    ax.axhline(y=0.4, color='orange', linestyle='--', alpha=0.5, linewidth=2)
    ax.set_xlabel('Mean Tanimoto', fontsize=11, fontweight='bold')
    ax.set_ylabel('Max Tanimoto', fontsize=11, fontweight='bold')
    ax.set_title('Mean vs Max Tanimoto', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax, label='Max Tanimoto')
    ax.grid(alpha=0.3)
    
    if has_score:
        # 5. Score vs AD
        ax = axes[4]
        for status in df['AD_Status'].unique():
            subset = df[df['AD_Status'] == status]
            ax.scatter(subset['MaxTanimoto'], subset['Total_Score'], 
                      label=status, alpha=0.6, s=30)
        ax.set_xlabel('Max Tanimoto', fontsize=11, fontweight='bold')
        ax.set_ylabel('Total Score', fontsize=11, fontweight='bold')
        ax.set_title('Score vs AD', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        
        # 6. Risk Score分布
        ax = axes[5]
        risk_counts = df.groupby(['AD_Status', 'Risk_Score']).size().unstack(fill_value=0)
        risk_counts.plot(kind='bar', stacked=True, ax=ax, 
                        color=['#2ecc71', '#f39c12', '#e74c3c'])
        ax.set_xlabel('AD Status', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title('Risk Distribution', fontsize=12, fontweight='bold')
        ax.legend(title='Risk Score', fontsize=8)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ 图表已保存: {save_path}")
    plt.close()


def generate_text_report(df, output_path):
    """生成文本报告"""
    
    total = len(df)
    inside = (df['AD_Status'] == 'Inside (Safe)').sum()
    boundary = (df['AD_Status'] == 'Boundary (Caution)').sum()
    outside = (df['AD_Status'] == 'Outside (High Risk)').sum()
    
    high_risk = df[df['MaxTanimoto'] < 0.4].sort_values('MaxTanimoto')
    
    has_score = 'Total_Score' in df.columns
    
    report = []
    report.append("="*80)
    report.append(f"{PROJECT_NAME}")
    report.append("适用域(AD)诊断报告")
    report.append("="*80)
    report.append("")
    report.append("## 1. 概览")
    report.append(f"   总分子数: {total}")
    report.append(f"   分析时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    report.append("## 2. AD分布")
    report.append(f"   ✅ Inside (Safe):         {inside:6d} ({inside/total*100:5.1f}%)")
    report.append(f"   ⚠️  Boundary (Caution):    {boundary:6d} ({boundary/total*100:5.1f}%)")
    report.append(f"   🚨 Outside (High Risk):   {outside:6d} ({outside/total*100:5.1f}%)")
    report.append("")
    
    report.append("## 3. 相似度统计")
    report.append(f"   平均 Max Tanimoto: {df['MaxTanimoto'].mean():.4f}")
    report.append(f"   中位数:            {df['MaxTanimoto'].median():.4f}")
    report.append(f"   最小值:            {df['MaxTanimoto'].min():.4f}")
    report.append(f"   最大值:            {df['MaxTanimoto'].max():.4f}")
    report.append(f"   标准差:            {df['MaxTanimoto'].std():.4f}")
    report.append("")
    
    if has_score:
        report.append("## 4. Score分析")
        report.append(f"   平均 Total Score:  {df['Total_Score'].mean():.4f}")
        report.append(f"   中位数:            {df['Total_Score'].median():.4f}")
        report.append(f"   最高分:            {df['Total_Score'].max():.4f}")
        report.append("")
        
        # 高分高风险分子
        high_score_high_risk = df[(df['Total_Score'] > 0.7) & (df['Risk_Score'] >= 3)]
        report.append(f"   ⚠️  高分但高风险分子: {len(high_score_high_risk)} 个")
        report.append("      （这些分子score高但可能不可靠，需谨慎对待）")
        report.append("")
    
    report.append("## 5. 风险评估")
    
    outside_pct = outside / total * 100
    if outside_pct > 30:
        report.append(f"   🚨 高风险警告: {outside_pct:.1f}% 分子在AD外！")
        report.append("")
        report.append("   ⚠️  建议:")
        report.append("   1. 这些分子的QSAR预测可能完全不准确")
        report.append("   2. 建议补充训练数据覆盖这些化学空间")
        report.append("   3. 对AD外的高分分子进行实验验证前要格外谨慎")
        report.append("   4. 考虑增加多样性惩罚到RL reward函数")
    elif outside_pct > 10:
        report.append(f"   ⚠️  中等风险: {outside_pct:.1f}% 分子在AD外")
        report.append("   建议对这些分子的预测保持谨慎态度")
    else:
        report.append(f"   ✅ 低风险: 仅{outside_pct:.1f}% 分子在AD外")
        report.append("   大部分分子的QSAR预测较为可靠")
    
    report.append("")
    report.append("## 6. 高风险分子详情")
    if len(high_risk) > 0:
        report.append(f"   共 {len(high_risk)} 个分子的 Max Tanimoto < 0.4")
        report.append("")
        report.append("   最高风险的10个分子:")
        for i, (idx, row) in enumerate(high_risk.head(10).iterrows(), 1):
            score_info = ""
            if has_score:
                score_info = f"  Score={row['Total_Score']:.3f}"
            report.append(f"   {i:2d}. {row['SMILES'][:55]:55s}  Sim={row['MaxTanimoto']:.4f}{score_info}")
    else:
        report.append("   ✅ 无高风险分子")
    
    report.append("")
    report.append("="*80)
    report.append("📊 详细结果文件:")
    report.append(f"  - {OUTPUT_DIR}/ad_analysis_results.csv         (完整AD分析)")
    report.append(f"  - {OUTPUT_DIR}/high_risk_molecules.csv         (高风险分子)")
    report.append(f"  - {OUTPUT_DIR}/ad_distribution.png             (可视化图表)")
    
    if has_score:
        report.append(f"  - {OUTPUT_DIR}/high_score_high_risk.csv        (高分高风险)")
    
    report.append("="*80)
    
    # 保存报告
    report_text = "\n".join(report)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print("\n" + report_text)


# ============================================================================
# 主函数
# ============================================================================

def main():
    """主流程"""
    
    print("="*80)
    print(f"{PROJECT_NAME}")
    print("适用域(AD)诊断")
    print("="*80)
    print()
    
    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Step 1: 读取训练集
    print("Step 1/4: 读取训练集")
    print("-" * 80)
    training_smiles = load_training_data()
    
    if training_smiles is None or len(training_smiles) == 0:
        print("✗ 无法加载训练集，终止")
        return
    
    print()
    
    # Step 2: 读取生成分子
    print("Step 2/4: 读取生成分子")
    print("-" * 80)
    generated_smiles, extra_cols = load_generated_data()
    
    if generated_smiles is None or len(generated_smiles) == 0:
        print("✗ 无法加载生成分子，终止")
        return
    
    print()
    
    # Step 3: 计算指纹和AD分析
    print("Step 3/4: 计算指纹和AD分析")
    print("-" * 80)
    
    train_fps, train_valid = compute_fingerprints(training_smiles)
    print(f"✓ 训练集有效分子: {len(train_valid)}/{len(training_smiles)}")
    
    gen_fps, gen_valid = compute_fingerprints(generated_smiles)
    print(f"✓ 生成集有效分子: {len(gen_valid)}/{len(generated_smiles)}")
    
    # 如果有额外列，需要对应筛选
    if extra_cols:
        filtered_extra = {}
        valid_indices = [i for i, smi in enumerate(generated_smiles) if Chem.MolFromSmiles(smi) is not None]
        for key, values in extra_cols.items():
            filtered_extra[key] = [values[i] for i in valid_indices]
        extra_cols = filtered_extra
    
    ad_results = calculate_ad_analysis(train_fps, gen_valid, gen_fps, extra_cols)
    
    # 保存结果
    result_path = os.path.join(OUTPUT_DIR, "ad_analysis_results.csv")
    ad_results.to_csv(result_path, index=False)
    print(f"\n✓ 完整结果已保存: {result_path}")
    
    # 高风险分子
    high_risk = ad_results[ad_results['MaxTanimoto'] < 0.4].sort_values('MaxTanimoto')
    high_risk_path = os.path.join(OUTPUT_DIR, "high_risk_molecules.csv")
    high_risk.to_csv(high_risk_path, index=False)
    print(f"✓ 高风险分子已保存: {high_risk_path} ({len(high_risk)} 个)")
    
    # 高分高风险分子
    if 'Total_Score' in ad_results.columns:
        high_score_high_risk = ad_results[
            (ad_results['Total_Score'] > 0.7) & (ad_results['Risk_Score'] >= 3)
        ].sort_values('Total_Score', ascending=False)
        
        hs_hr_path = os.path.join(OUTPUT_DIR, "high_score_high_risk.csv")
        high_score_high_risk.to_csv(hs_hr_path, index=False)
        print(f"✓ 高分高风险分子已保存: {hs_hr_path} ({len(high_score_high_risk)} 个)")
    
    print()
    
    # Step 4: 可视化和报告
    print("Step 4/4: 生成可视化和报告")
    print("-" * 80)
    
    plot_path = os.path.join(OUTPUT_DIR, "ad_distribution.png")
    plot_ad_distribution(ad_results, plot_path)
    
    report_path = os.path.join(OUTPUT_DIR, "diagnostic_report.txt")
    generate_text_report(ad_results, report_path)
    print(f"✓ 报告已保存: {report_path}")
    
    print()
    print("="*80)
    print("✅✅✅ 诊断完成！✅✅✅")
    print("="*80)
    print(f"\n📁 所有结果保存在: {OUTPUT_DIR}/")
    print("\n🔍 关键发现:")
    
    inside = (ad_results['AD_Status'] == 'Inside (Safe)').sum()
    outside = (ad_results['AD_Status'] == 'Outside (High Risk)').sum()
    total = len(ad_results)
    
    print(f"   • {inside}/{total} ({inside/total*100:.1f}%) 分子在AD内（可靠）")
    print(f"   • {outside}/{total} ({outside/total*100:.1f}%) 分子在AD外（高风险）")
    
    if 'Total_Score' in ad_results.columns:
        high_score_count = (ad_results['Total_Score'] > 0.7).sum()
        print(f"   • {high_score_count} 个高分分子 (Score > 0.7)")
    
    print()


if __name__ == "__main__":
    main()