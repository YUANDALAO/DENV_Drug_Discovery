#!/usr/bin/env python3
"""
评估AlphaFold预测的CIF文件质量
并与4M9K实验结构比对
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    from Bio.PDB import MMCIFParser, PDBParser, Superimposer, PDBIO
except ImportError:
    print("请安装Biopython: pip install biopython")
    exit(1)


class AlphaFoldCIFEvaluator:
    """AlphaFold CIF文件评估器"""
    
    def __init__(self, cif_file: str):
        self.cif_file = cif_file
        self.parser = MMCIFParser(QUIET=True)
        self.structure = None
        self.plddt_scores = []
        self.residue_numbers = []
        
    def load_structure(self):
        """加载CIF结构"""
        print(f"\n📂 加载AlphaFold预测结构...")
        print(f"   文件: {Path(self.cif_file).name}")
        
        try:
            self.structure = self.parser.get_structure('AF', self.cif_file)
            print("   ✓ 结构加载成功")
            return True
        except Exception as e:
            print(f"   ✗ 加载失败: {e}")
            return False
    
    def extract_plddt(self):
        """从CIF文件提取pLDDT分数"""
        print("\n📊 提取pLDDT分数...")
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    try:
                        ca = residue['CA']
                        plddt = ca.bfactor
                        res_num = residue.id[1]
                        
                        self.plddt_scores.append(plddt)
                        self.residue_numbers.append(res_num)
                    except:
                        continue
        
        if not self.plddt_scores:
            print("   ✗ 未找到pLDDT分数")
            return False
        
        self.plddt_scores = np.array(self.plddt_scores)
        self.residue_numbers = np.array(self.residue_numbers)
        
        print(f"   ✓ 成功提取 {len(self.plddt_scores)} 个残基的pLDDT")
        return True
    
    def analyze_quality(self):
        """分析预测质量"""
        print("\n" + "="*70)
        print("AlphaFold 预测质量分析")
        print("="*70)
        
        mean_plddt = self.plddt_scores.mean()
        min_plddt = self.plddt_scores.min()
        max_plddt = self.plddt_scores.max()
        
        print(f"\n整体质量:")
        print(f"  平均 pLDDT: {mean_plddt:.2f}")
        print(f"  最低 pLDDT: {min_plddt:.2f}")
        print(f"  最高 pLDDT: {max_plddt:.2f}")
        
        # 分布统计
        excellent = (self.plddt_scores > 90).sum()
        good = ((self.plddt_scores > 70) & (self.plddt_scores <= 90)).sum()
        moderate = ((self.plddt_scores > 50) & (self.plddt_scores <= 70)).sum()
        poor = (self.plddt_scores <= 50).sum()
        total = len(self.plddt_scores)
        
        print(f"\npLDDT分布:")
        print(f"  >90 (优秀):   {excellent:3d}/{total} ({excellent/total*100:5.1f}%)")
        print(f"  70-90 (良好): {good:3d}/{total} ({good/total*100:5.1f}%)")
        print(f"  50-70 (一般): {moderate:3d}/{total} ({moderate/total*100:5.1f}%)")
        print(f"  <50 (差):     {poor:3d}/{total} ({poor/total*100:5.1f}%)")
        
        # 活性位点区域分析
        active_site_mask = (self.residue_numbers >= 50) & (self.residue_numbers <= 150)
        if active_site_mask.any():
            active_plddt = self.plddt_scores[active_site_mask]
            print(f"\n活性位点区域 (残基 50-150):")
            print(f"  平均 pLDDT: {active_plddt.mean():.2f}")
            print(f"  最低 pLDDT: {active_plddt.min():.2f}")
            
            if active_plddt.mean() > 80:
                print("  ✓ 活性位点质量优秀")
                quality = "excellent"
            elif active_plddt.mean() > 70:
                print("  ⚠ 活性位点质量良好")
                quality = "good"
            else:
                print("  ✗ 活性位点质量不佳")
                quality = "poor"
        else:
            quality = "unknown"
        
        # 整体评价
        print("\n" + "-"*70)
        if mean_plddt > 85:
            print("✓ 整体质量: 优秀 - 可用于高精度对接")
            overall = "excellent"
        elif mean_plddt > 75:
            print("⚠ 整体质量: 良好 - 建议与实验结构对比验证")
            overall = "good"
        elif mean_plddt > 65:
            print("⚠ 整体质量: 一般 - 建议优先使用实验结构")
            overall = "moderate"
        else:
            print("✗ 整体质量: 较差 - 不建议用于对接")
            overall = "poor"
        print("="*70)
        
        return {
            'mean': mean_plddt,
            'min': min_plddt,
            'max': max_plddt,
            'overall_quality': overall,
            'active_site_quality': quality
        }
    
    def plot_plddt(self, output_file='plddt_analysis.png'):
        """绘制pLDDT分布图"""
        print(f"\n📊 生成pLDDT分布图...")
        
        # 设置字体
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
        plt.rcParams['axes.unicode_minus'] = False
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
        
        # 上图：沿序列的pLDDT
        ax1.plot(self.residue_numbers, self.plddt_scores, linewidth=1.5, color='#2563eb')
        ax1.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='Excellent (>90)')
        ax1.axhline(y=70, color='orange', linestyle='--', alpha=0.5, label='Good (>70)')
        ax1.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='Acceptable (>50)')
        ax1.fill_between(self.residue_numbers, self.plddt_scores, 50, 
                         where=(self.plddt_scores >= 50), alpha=0.3, color='lightgreen')
        ax1.fill_between(self.residue_numbers, self.plddt_scores, 50, 
                         where=(self.plddt_scores < 50), alpha=0.3, color='lightcoral')
        
        ax1.set_xlabel('Residue Number', fontsize=12)
        ax1.set_ylabel('pLDDT Score', fontsize=12)
        ax1.set_title(f'AlphaFold Prediction Quality (Mean pLDDT: {self.plddt_scores.mean():.1f})', 
                     fontsize=14, fontweight='bold')
        ax1.legend(loc='lower left')
        ax1.grid(alpha=0.3)
        ax1.set_ylim([0, 100])
        
        # 标注活性位点区域
        ax1.axvspan(50, 150, alpha=0.1, color='purple')
        ax1.text(100, 5, 'Active Site Region', ha='center', fontsize=10, 
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # 下图：pLDDT直方图
        ax2.hist(self.plddt_scores, bins=50, color='#3b82f6', edgecolor='black', alpha=0.7)
        ax2.axvline(x=90, color='green', linestyle='--', linewidth=2, label='Excellent threshold')
        ax2.axvline(x=70, color='orange', linestyle='--', linewidth=2, label='Good threshold')
        ax2.axvline(x=self.plddt_scores.mean(), color='red', linestyle='-', 
                   linewidth=2, label=f'Mean ({self.plddt_scores.mean():.1f})')
        
        ax2.set_xlabel('pLDDT Score', fontsize=12)
        ax2.set_ylabel('Residue Count', fontsize=12)
        ax2.set_title('pLDDT Distribution', fontsize=14, fontweight='bold')
        ax2.legend()
        ax2.grid(alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"   ✓ 图表已保存: {output_file}")
        
        return output_file
    
    def convert_to_pdb(self, output_file='alphafold_converted.pdb'):
        """将CIF转换为PDB格式"""
        print(f"\n💾 转换为PDB格式...")
        
        try:
            io = PDBIO()
            io.set_structure(self.structure)
            io.save(output_file)
            print(f"   ✓ PDB已保存: {output_file}")
            return output_file
        except Exception as e:
            print(f"   ✗ 转换失败: {e}")
            return None


def compare_with_experimental(af_pdb, exp_pdb='structures/4M9K.pdb'):
    """与实验结构比对 - 完全修复版本"""
    print("\n" + "="*70)
    print("与实验结构(4M9K)比对")
    print("="*70)
    
    parser = PDBParser(QUIET=True)
    
    try:
        af_struct = parser.get_structure('AF', af_pdb)
        exp_struct = parser.get_structure('4M9K', exp_pdb)
    except Exception as e:
        print(f"✗ 加载结构失败: {e}")
        return None
    
    # 先检查残基范围并收集所有存在的残基
    print("\n检查残基编号范围...")
    
    af_residues = {}
    for res in af_struct[0]['A']:
        if res.id[0] == ' ':  # 只要标准残基
            try:
                if 'CA' in res:
                    af_residues[res.id[1]] = res
            except:
                continue
    
    exp_residues = {}
    for res in exp_struct[0]['A']:
        if res.id[0] == ' ':
            try:
                if 'CA' in res:
                    exp_residues[res.id[1]] = res
            except:
                continue
    
    af_res_nums = sorted(af_residues.keys())
    exp_res_nums = sorted(exp_residues.keys())
    
    print(f"  AlphaFold: {af_res_nums[0]} - {af_res_nums[-1]} (共{len(af_res_nums)}个)")
    print(f"  4M9K: {exp_res_nums[0]} - {exp_res_nums[-1]} (共{len(exp_res_nums)}个)")
    
    # 策略：从AlphaFold的第1个残基开始，对应4M9K的第62个残基
    # 因为AlphaFold预测的是完整NS3蛋白酶域(1-185)
    # 而4M9K包含NS2B(1-47) + linker(48-61) + NS3(62-245)
    
    af_start = 1
    af_end = min(185, len(af_res_nums))  # NS3蛋白酶域最多185个残基
    
    exp_start = 62  # 4M9K中NS3的起始
    
    print(f"\n对齐策略:")
    print(f"  AlphaFold残基 {af_start} → 4M9K残基 {exp_start}")
    print(f"  (AlphaFold是完整NS3，4M9K的NS3从62开始)")
    
    # 提取匹配的Cα原子对
    af_atoms = []
    exp_atoms = []
    matched_pairs = []
    
    for i in range(af_end):
        af_res_num = af_start + i
        exp_res_num = exp_start + i
        
        # 检查两个残基是否都存在
        if af_res_num in af_residues and exp_res_num in exp_residues:
            try:
                af_ca = af_residues[af_res_num]['CA']
                exp_ca = exp_residues[exp_res_num]['CA']
                
                af_atoms.append(af_ca)
                exp_atoms.append(exp_ca)
                matched_pairs.append((af_res_num, exp_res_num))
            except:
                continue
    
    num_matched = len(af_atoms)
    
    if num_matched < 50:
        print(f"\n⚠ 警告: 可比对的残基太少 ({num_matched})")
        print("   无法进行可靠的RMSD计算")
        print("\n   可能原因:")
        print("   - AlphaFold预测的序列与4M9K不完全对应")
        print("   - 结构中有缺失残基")
        return None
    
    print(f"\n✓ 成功匹配 {num_matched} 个残基对")
    print(f"  对应范围: AF {matched_pairs[0][0]}-{matched_pairs[-1][0]} "
          f"↔ 4M9K {matched_pairs[0][1]}-{matched_pairs[-1][1]}")
    
    # 执行叠合
    print(f"\n正在叠合 {num_matched} 个Cα原子...")
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(exp_atoms, af_atoms)
    super_imposer.apply(af_struct.get_atoms())
    
    rmsd = super_imposer.rms
    
    print(f"\n结果:")
    print(f"  整体RMSD: {rmsd:.2f} Å (基于{num_matched}个残基)")
    
    # 评估
    if rmsd < 2.0:
        print("  ✓ 结构高度相似 - AlphaFold预测优秀")
        quality = "excellent"
    elif rmsd < 3.0:
        print("  ⚠ 结构相似度良好 - 可以使用")
        quality = "good"
    elif rmsd < 4.0:
        print("  ⚠ 结构有一定差异 - 建议双轨验证")
        quality = "moderate"
    else:
        print("  ✗ 结构差异较大 - 建议优先用实验结构")
        quality = "poor"
    
    # 保存叠合后的结构
    aligned_file = 'structures/alphafold_aligned_to_4M9K.pdb'
    io = PDBIO()
    io.set_structure(af_struct)
    io.save(aligned_file)
    print(f"\n✓ 叠合结构已保存: {aligned_file}")
    print("  可用PyMOL查看:")
    print("    pymol structures/4M9K.pdb structures/alphafold_aligned_to_4M9K.pdb")
    
    return {
        'rmsd': rmsd,
        'quality': quality,
        'aligned_file': aligned_file,
        'num_residues': num_matched
    }


def generate_recommendation(quality_results, comparison_results):
    """生成使用建议"""
    print("\n" + "="*70)
    print("🎯 使用建议")
    print("="*70)
    
    af_quality = quality_results['overall_quality']
    active_quality = quality_results['active_site_quality']
    
    if comparison_results:
        rmsd_quality = comparison_results['quality']
        rmsd = comparison_results['rmsd']
        num_res = comparison_results.get('num_residues', 0)
    else:
        rmsd_quality = 'unknown'
        rmsd = None
        num_res = 0
    
    print(f"\n评估结果:")
    print(f"  AlphaFold整体质量: {af_quality}")
    print(f"  活性位点质量: {active_quality}")
    if rmsd is not None:
        print(f"  与4M9K相似度: {rmsd_quality} (RMSD={rmsd:.2f}Å, {num_res}残基)")
    
    print("\n" + "-"*70)
    
    # 决策逻辑
    if af_quality == 'excellent' and active_quality == 'excellent' and rmsd_quality in ['excellent', 'good']:
        print("✅ 推荐方案: 优先使用AlphaFold模型")
        print("\n理由:")
        print("  • 预测质量优秀(pLDDT>85)")
        print("  • 活性位点可靠(pLDDT>80)")
        print("  • 与实验结构高度一致(RMSD<3Å)")
        print("  • 序列完全匹配你的训练数据(P29991)")
        print("\nREINVENT4配置:")
        print("  QSAR: 0.6")
        print("  Docking(AlphaFold): 0.3")
        print("  QED: 0.1")
        
    elif (af_quality in ['good', 'excellent'] and active_quality in ['good', 'excellent'] 
          and rmsd_quality in ['excellent', 'good', 'moderate']):
        print("⚠️ 推荐方案: 双轨验证")
        print("\n理由:")
        print("  • AlphaFold质量良好且活性位点优秀")
        print("  • 与实验结构相似度可接受")
        print("  • 同时使用两个结构更稳妥")
        print("\nREINVENT4配置:")
        print("  QSAR: 0.5")
        print("  Docking(4M9K): 0.25")
        print("  Docking(AlphaFold): 0.15")
        print("  QED: 0.1")
        
    else:
        print("🔴 推荐方案: 优先使用4M9K实验结构")
        print("\n理由:")
        print("  • AlphaFold预测质量不够理想")
        print("  • 或与实验结构差异较大")
        print("  • 4M9K是高分辨率实验结构(1.46Å)")
        print("\nREINVENT4配置:")
        print("  QSAR: 0.6")
        print("  Docking(4M9K): 0.3")
        print("  QED: 0.1")
    
    print("="*70)


def main():
    """主函数"""
    import sys
    
    # 设置CIF文件路径
    if len(sys.argv) < 2:
        # 默认路径（WSL格式）
        cif_file = "/mnt/c/Users/ucsaheu/python_projects/DENV_Drug_Discovery/05_Generative_AI_REINVENT/REINVENT4-main/denv_inhibitor_design/scripts/fold_2025_10_15_16_12_model_0.cif"
    else:
        cif_file = sys.argv[1]
    
    # 检查文件是否存在
    if not Path(cif_file).exists():
        print(f"\n✗ 错误: 找不到文件 {cif_file}")
        print("\n使用方法:")
        print("  python evaluate_alphafold_cif.py <cif文件路径>")
        return
    
    print("\n" + "="*70)
    print("AlphaFold CIF 文件评估工具")
    print("="*70)
    
    # 创建必要的目录
    Path("results").mkdir(exist_ok=True)
    Path("structures").mkdir(exist_ok=True)
    
    # 1. 加载和评估AlphaFold结构
    evaluator = AlphaFoldCIFEvaluator(cif_file)
    
    if not evaluator.load_structure():
        return
    
    if not evaluator.extract_plddt():
        return
    
    quality_results = evaluator.analyze_quality()
    
    # 2. 绘制pLDDT图
    try:
        evaluator.plot_plddt('results/alphafold_plddt_analysis.png')
    except Exception as e:
        print(f"⚠ 绘图失败: {e}")
        print("  继续执行其他步骤...")
    
    # 3. 转换为PDB
    af_pdb = evaluator.convert_to_pdb('structures/DENV2_alphafold.pdb')
    
    # 4. 与实验结构比对（如果4M9K存在）
    exp_pdb = 'structures/4M9K.pdb'
    if Path(exp_pdb).exists() and af_pdb:
        comparison_results = compare_with_experimental(af_pdb, exp_pdb)
    else:
        print(f"\n⚠ 未找到4M9K.pdb，跳过结构比对")
        print(f"   建议下载: wget https://files.rcsb.org/download/4M9K.pdb -O {exp_pdb}")
        comparison_results = None
    
    # 5. 生成建议
    generate_recommendation(quality_results, comparison_results)
    
    print("\n✓ 评估完成!")
    print(f"   - pLDDT图表: results/alphafold_plddt_analysis.png")
    print(f"   - PDB文件: structures/DENV2_alphafold.pdb")
    if comparison_results:
        print(f"   - 叠合结构: {comparison_results['aligned_file']}")
    
    print("\n下一步:")
    if not Path(exp_pdb).exists():
        print("  1. 下载4M9K: cd structures && wget https://files.rcsb.org/download/4M9K.pdb")
        print("  2. 重新运行脚本查看RMSD比对结果")
    else:
        print("  1. 根据上述建议配置REINVENT4")
        print("  2. 准备对接配置文件")
        print("  3. 开始生成式设计！")


if __name__ == "__main__":
    main()