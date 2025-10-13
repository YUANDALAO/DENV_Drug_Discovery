#!/usr/bin/env python3
"""
统一的候选物分析和可视化工具
结合筛选、排序、统计和结构可视化

用法: python unified_candidate_analysis.py experiments/runs/run9_t1200
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import seaborn as sns

# 配置matplotlib中文显示
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

class CandidateAnalyzer:
    """候选物分析器"""
    
    def __init__(self, results_csv, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"📂 读取数据: {results_csv}")
        self.df = pd.read_csv(results_csv)
        
        # 自动检测列名
        self._detect_columns()
        
        # 计算IC50
        if self.activity_col:
            self.df['IC50_nM'] = 10 ** (9 - self.df[self.activity_col])
        
        # 去重
        self.df_unique = self.df.drop_duplicates(subset=['SMILES'])
        print(f"✓ 总分子数: {len(self.df):,}, 去重后: {len(self.df_unique):,}")
    
    def _detect_columns(self):
        """自动检测列名"""
        def find_col(patterns):
            for pattern in patterns:
                matches = [col for col in self.df.columns if pattern.lower() in col.lower()]
                if matches:
                    return matches[0]
            return None
        
        self.activity_col = find_col(['DENV_Activity (raw)', 'Activity (raw)', 'pIC50'])
        self.qed_col = find_col(['QED (raw)', 'Drug_Likeness', 'QED'])
        self.sa_col = find_col(['SA (raw)', 'Synthetic_Accessibility', 'SA'])
        self.mw_col = find_col(['MW (raw)', 'Molecular weight', 'MW'])
        self.logp_col = find_col(['LogP (raw)', 'SlogP', 'LogP'])
        
        print(f"\n📊 检测到的列:")
        print(f"  • Activity: {self.activity_col}")
        print(f"  • QED: {self.qed_col}")
        print(f"  • SA: {self.sa_col}")
        print(f"  • MW: {self.mw_col}")
        print(f"  • LogP: {self.logp_col}")
    
    def filter_by_criteria(self, criteria_name, min_pic50, min_qed, max_sa, 
                          mw_range, logp_range):
        """按条件筛选"""
        mask = pd.Series(True, index=self.df_unique.index)
        
        if self.activity_col:
            mask &= self.df_unique[self.activity_col] >= min_pic50
        if self.qed_col:
            mask &= self.df_unique[self.qed_col] >= min_qed
        if self.sa_col:
            mask &= self.df_unique[self.sa_col] <= max_sa
        if self.mw_col:
            mask &= (self.df_unique[self.mw_col] >= mw_range[0]) & \
                    (self.df_unique[self.mw_col] <= mw_range[1])
        if self.logp_col:
            mask &= (self.df_unique[self.logp_col] >= logp_range[0]) & \
                    (self.df_unique[self.logp_col] <= logp_range[1])
        
        filtered = self.df_unique[mask].copy()
        filtered = filtered.sort_values('Score', ascending=False)
        
        return filtered
    
    def analyze_multi_criteria(self):
        """多层次筛选分析"""
        print("\n" + "="*70)
        print("🔍 多层次候选物筛选".center(70))
        print("="*70)
        
        criteria_sets = [
            {
                'name': '💎 金标准 (Gold Standard)',
                'short': 'gold',
                'min_pic50': 8.0,
                'min_qed': 0.7,
                'max_sa': 4.0,
                'mw_range': (300, 500),
                'logp_range': (1, 4)
            },
            {
                'name': '⭐ 高标准 (High Quality)',
                'short': 'high',
                'min_pic50': 7.5,
                'min_qed': 0.6,
                'max_sa': 4.5,
                'mw_range': (250, 550),
                'logp_range': (0.5, 5)
            },
            {
                'name': '✓ 中标准 (Good Quality)',
                'short': 'good',
                'min_pic50': 7.0,
                'min_qed': 0.5,
                'max_sa': 5.0,
                'mw_range': (200, 600),
                'logp_range': (0, 6)
            }
        ]
        
        results = {}
        
        for criteria in criteria_sets:
            print(f"\n{'='*70}")
            print(f"{criteria['name']}")
            print(f"{'='*70}")
            print(f"筛选条件: pIC50≥{criteria['min_pic50']}, QED≥{criteria['min_qed']}, "
                  f"SA≤{criteria['max_sa']}")
            print(f"           MW: {criteria['mw_range']}, LogP: {criteria['logp_range']}")
            
            candidates = self.filter_by_criteria(
                criteria['name'],
                criteria['min_pic50'],
                criteria['min_qed'],
                criteria['max_sa'],
                criteria['mw_range'],
                criteria['logp_range']
            )
            
            results[criteria['short']] = candidates
            
            print(f"\n✓ 符合条件: {len(candidates)} 个分子 "
                  f"({len(candidates)/len(self.df_unique)*100:.2f}%)")
            
            if len(candidates) > 0:
                self._print_statistics(candidates)
                
                # 保存CSV
                output_file = self.output_dir / f"candidates_{criteria['short']}.csv"
                candidates.to_csv(output_file, index=False)
                print(f"\n💾 已保存: {output_file}")
                
                # 显示Top 5
                print(f"\n🏆 Top 5 分子:")
                self._print_top_molecules(candidates, n=5)
        
        return results
    
    def _print_statistics(self, df):
        """打印统计信息"""
        print(f"\n📈 统计范围:")
        
        if self.activity_col and self.activity_col in df.columns:
            pic50_vals = df[self.activity_col].dropna()
            if len(pic50_vals) > 0:
                print(f"  • pIC50: {pic50_vals.min():.2f} - {pic50_vals.max():.2f} "
                      f"(平均: {pic50_vals.mean():.2f})")
                ic50_vals = 10**(9 - pic50_vals)
                print(f"  • IC50: {ic50_vals.min():.1f} - {ic50_vals.max():.1f} nM "
                      f"(平均: {ic50_vals.mean():.1f} nM)")
        
        if self.qed_col and self.qed_col in df.columns:
            qed_vals = df[self.qed_col].dropna()
            if len(qed_vals) > 0:
                print(f"  • QED: {qed_vals.min():.2f} - {qed_vals.max():.2f} "
                      f"(平均: {qed_vals.mean():.2f})")
        
        if self.sa_col and self.sa_col in df.columns:
            sa_vals = df[self.sa_col].dropna()
            if len(sa_vals) > 0:
                print(f"  • SA: {sa_vals.min():.2f} - {sa_vals.max():.2f} "
                      f"(平均: {sa_vals.mean():.2f})")
        
        print(f"  • 总分: {df['Score'].min():.3f} - {df['Score'].max():.3f} "
              f"(平均: {df['Score'].mean():.3f})")
    
    def _print_top_molecules(self, df, n=5):
        """打印Top分子"""
        for idx, (i, row) in enumerate(df.head(n).iterrows(), 1):
            print(f"\n  {idx}. Score: {row['Score']:.3f}")
            print(f"     SMILES: {row['SMILES']}")
            
            details = []
            if self.activity_col and self.activity_col in row:
                pic50 = row[self.activity_col]
                ic50 = 10**(9 - pic50)
                details.append(f"pIC50={pic50:.2f} (IC50={ic50:.1f}nM)")
            if self.qed_col and self.qed_col in row:
                details.append(f"QED={row[self.qed_col]:.2f}")
            if self.sa_col and self.sa_col in row:
                details.append(f"SA={row[self.sa_col]:.2f}")
            
            if details:
                print(f"     {', '.join(details)}")
    
    def visualize_candidates(self, candidates_dict, n_per_criteria=10):
        """可视化各标准的Top候选物"""
        print("\n" + "="*70)
        print("🎨 生成结构可视化".center(70))
        print("="*70)
        
        for criteria_name, df in candidates_dict.items():
            if len(df) == 0:
                print(f"\n⚠️  {criteria_name}: 无候选物，跳过可视化")
                continue
            
            n_mols = min(n_per_criteria, len(df))
            print(f"\n📸 {criteria_name}: 生成 Top {n_mols} 分子结构图...")
            
            top_df = df.head(n_mols)
            
            # 解析分子
            mols = []
            legends = []
            
            for _, row in top_df.iterrows():
                mol = Chem.MolFromSmiles(row['SMILES'])
                if mol:
                    mols.append(mol)
                    
                    # 构建图例
                    legend_parts = [f"Score: {row['Score']:.3f}"]
                    
                    if self.activity_col and self.activity_col in row:
                        pic50 = row[self.activity_col]
                        ic50 = 10**(9 - pic50)
                        legend_parts.append(f"pIC50: {pic50:.2f}")
                        legend_parts.append(f"IC50: {ic50:.1f} nM")
                    
                    if self.qed_col and self.qed_col in row:
                        legend_parts.append(f"QED: {row[self.qed_col]:.2f}")
                    
                    if self.sa_col and self.sa_col in row:
                        legend_parts.append(f"SA: {row[self.sa_col]:.2f}")
                    
                    legends.append("\n".join(legend_parts))
            
            if not mols:
                print(f"  ⚠️  无有效分子结构")
                continue
            
            # 生成图像
            n_cols = min(4, len(mols))
            n_rows = (len(mols) + n_cols - 1) // n_cols
            
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=n_cols,
                subImgSize=(400, 400),
                legends=legends,
                returnPNG=False
            )
            
            output_file = self.output_dir / f"structures_{criteria_name}_top{n_mols}.png"
            img.save(output_file)
            print(f"  ✓ 已保存: {output_file}")
    
    def compare_scoring_methods(self):
        """对比不同排序方法的结果"""
        print("\n" + "="*70)
        print("📊 对比分析: 总分排序 vs 多维筛选".center(70))
        print("="*70)
        
        # 方法1: 按总分排序Top 20
        top_by_score = self.df_unique.nlargest(20, 'Score')
        
        # 方法2: 金标准筛选
        gold_candidates = self.filter_by_criteria(
            'gold', min_pic50=8.0, min_qed=0.7, max_sa=4.0,
            mw_range=(300, 500), logp_range=(1, 4)
        ).head(20)
        
        print("\n🔵 方法1: 按总分排序 (Top 20)")
        print("   - 优点: 综合平衡性好")
        print("   - 缺点: 可能某些关键属性不够优秀")
        self._print_statistics(top_by_score)
        
        if len(gold_candidates) > 0:
            print("\n🟡 方法2: 金标准多维筛选 (Top 20)")
            print("   - 优点: 每个维度都必须达标")
            print("   - 缺点: 筛选更严格，候选物更少")
            self._print_statistics(gold_candidates)
        
        # 交集分析
        common = set(top_by_score['SMILES']) & set(gold_candidates['SMILES'])
        print(f"\n🔄 两种方法的重叠: {len(common)} 个分子")
        
        if len(common) < 5:
            print("\n💡 建议: 总分Top分子的活性可能不够高！")
            print("   → 使用多维筛选方法可以找到活性更强的候选物")
    
    def generate_summary_report(self):
        """生成总结报告"""
        print("\n" + "="*70)
        print("📋 综合分析报告".center(70))
        print("="*70)
        
        # IC50分布
        if 'IC50_nM' in self.df_unique.columns:
            print("\n🎯 IC50活性分布:")
            
            ic50_ranges = [
                (0, 10, '0-10 nM (极高活性) 🔥'),
                (10, 50, '10-50 nM (高活性) ⭐'),
                (50, 100, '50-100 nM (中等活性) ✓'),
                (100, 1000, '100-1000 nM (低活性)'),
                (1000, float('inf'), '>1000 nM (极低活性)')
            ]
            
            for low, high, label in ic50_ranges:
                if high == float('inf'):
                    count = len(self.df_unique[self.df_unique['IC50_nM'] >= low])
                else:
                    count = len(self.df_unique[
                        (self.df_unique['IC50_nM'] >= low) & 
                        (self.df_unique['IC50_nM'] < high)
                    ])
                pct = count / len(self.df_unique) * 100
                bar = '█' * int(pct / 2)
                print(f"  {label:30s}: {count:5d} ({pct:5.2f}%) {bar}")
        
        # 多维度优秀分子
        print("\n🌟 多维度优秀分子:")
        
        if self.activity_col and self.qed_col and self.sa_col:
            excellent = self.df_unique[
                (self.df_unique[self.activity_col] > 7.5) |
                (self.df_unique[self.qed_col] > 0.75) |
                (self.df_unique[self.sa_col] < 3.5)
            ]
            print(f"  • 至少一个维度优秀: {len(excellent)} "
                  f"({len(excellent)/len(self.df_unique)*100:.1f}%)")
            
            excellent_2 = self.df_unique[
                ((self.df_unique[self.activity_col] > 7.5) & 
                 (self.df_unique[self.qed_col] > 0.65)) |
                ((self.df_unique[self.activity_col] > 7.5) & 
                 (self.df_unique[self.sa_col] < 4.0)) |
                ((self.df_unique[self.qed_col] > 0.65) & 
                 (self.df_unique[self.sa_col] < 4.0))
            ]
            print(f"  • 至少两个维度优秀: {len(excellent_2)} "
                  f"({len(excellent_2)/len(self.df_unique)*100:.1f}%)")
            
            excellent_3 = self.df_unique[
                (self.df_unique[self.activity_col] > 7.5) &
                (self.df_unique[self.qed_col] > 0.65) &
                (self.df_unique[self.sa_col] < 4.0)
            ]
            print(f"  • 三个维度都优秀: {len(excellent_3)} "
                  f"({len(excellent_3)/len(self.df_unique)*100:.1f}%)")
    
    def run_full_analysis(self):
        """运行完整分析流程"""
        # 1. 多标准筛选
        candidates_dict = self.analyze_multi_criteria()
        
        # 2. 可视化
        self.visualize_candidates(candidates_dict, n_per_criteria=10)
        
        # 3. 对比分析
        self.compare_scoring_methods()
        
        # 4. 总结报告
        self.generate_summary_report()
        
        print("\n" + "="*70)
        print("✅ 分析完成！".center(70))
        print("="*70)
        print(f"\n📁 所有结果已保存到: {self.output_dir}")


def main():
    if len(sys.argv) < 2:
        print("用法: python unified_candidate_analysis.py <实验文件夹>")
        print("示例: python unified_candidate_analysis.py experiments/runs/run9_t1200")
        sys.exit(1)
    
    run_dir = Path(sys.argv[1])
    
    # 查找结果文件
    results_files = list(run_dir.glob("results_*.csv"))
    
    if not results_files:
        print(f"❌ 错误: 在 {run_dir} 中找不到 results_*.csv 文件")
        sys.exit(1)
    
    # 使用最新的结果文件
    results_file = max(results_files, key=lambda x: x.stat().st_mtime)
    
    print("="*70)
    print("🔬 REINVENT4 候选物统一分析工具".center(70))
    print("="*70)
    print(f"\n📂 实验目录: {run_dir}")
    print(f"📄 结果文件: {results_file.name}\n")
    
    # 创建分析器并运行
    analyzer = CandidateAnalyzer(results_file, run_dir)
    analyzer.run_full_analysis()


if __name__ == "__main__":
    main()