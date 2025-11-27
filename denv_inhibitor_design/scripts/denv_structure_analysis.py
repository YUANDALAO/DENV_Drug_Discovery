#!/usr/bin/env python3
"""
DENV2 NS2B-NS3 Structure Prediction and Comparison with DENV3
自动化流程：AlphaFold预测 → 结构叠合 → 活性位点分析 → 决策建议
"""

import os
import requests
import json
from pathlib import Path
from Bio import SeqIO, Align
from Bio.PDB import PDBParser, Superimposer, PDBIO, Selection
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple

# ============================================================================
# 第1步：获取DENV2和DENV3的序列
# ============================================================================

class DENVSequenceFetcher:
    """从UniProt获取DENV序列"""
    
    UNIPROT_IDS = {
        'DENV2': 'P29991',  # DENV2 NS3 protease
        'DENV3': 'P27915',  # DENV3 NS3 protease
    }
    
    @staticmethod
    def fetch_sequence(uniprot_id: str, output_file: str) -> str:
        """从UniProt下载序列"""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url)
        
        if response.status_code == 200:
            with open(output_file, 'w') as f:
                f.write(response.text)
            print(f"✓ 已下载 {uniprot_id} 序列到 {output_file}")
            
            # 解析序列
            record = SeqIO.read(output_file, "fasta")
            return str(record.seq)
        else:
            raise Exception(f"✗ 无法下载 {uniprot_id}: {response.status_code}")
    
    @classmethod
    def fetch_all(cls, output_dir: str = "sequences"):
        """下载所有DENV序列"""
        Path(output_dir).mkdir(exist_ok=True)
        sequences = {}
        
        for serotype, uniprot_id in cls.UNIPROT_IDS.items():
            output_file = f"{output_dir}/{serotype}_NS3.fasta"
            sequences[serotype] = cls.fetch_sequence(uniprot_id, output_file)
        
        return sequences


# ============================================================================
# 第2步：序列比对分析
# ============================================================================

class SequenceAligner:
    """序列比对和保守性分析"""
    
    # DENV NS3蛋白酶催化三联体和活性位点关键残基（基于文献）
    CATALYTIC_TRIAD = [51, 75, 135]  # His51, Asp75, Ser135
    SUBSTRATE_BINDING = [36, 37, 51, 52, 82, 84, 85, 129, 130, 132, 133, 135, 152, 155]
    
    @staticmethod
    def align_sequences(seq1: str, seq2: str, label1: str, label2: str) -> Dict:
        """全局序列比对"""
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        
        identity = sum(a == b for a, b in zip(best_alignment[0], best_alignment[1]) 
                      if a != '-' and b != '-')
        length = min(len(seq1), len(seq2))
        identity_pct = (identity / length) * 100
        
        return {
            'alignment': best_alignment,
            'identity': identity,
            'length': length,
            'identity_percentage': identity_pct,
            'label1': label1,
            'label2': label2
        }
    
    @classmethod
    def analyze_active_site(cls, seq1: str, seq2: str) -> Dict:
        """分析活性位点残基保守性"""
        results = {
            'catalytic_triad': [],
            'binding_site': [],
            'triad_conserved': True,
            'binding_conserved_count': 0
        }
        
        # 检查催化三联体
        for pos in cls.CATALYTIC_TRIAD:
            if pos < len(seq1) and pos < len(seq2):
                res1, res2 = seq1[pos], seq2[pos]
                conserved = (res1 == res2)
                results['catalytic_triad'].append({
                    'position': pos,
                    'DENV2': res1,
                    'DENV3': res2,
                    'conserved': conserved
                })
                if not conserved:
                    results['triad_conserved'] = False
        
        # 检查底物结合位点
        for pos in cls.SUBSTRATE_BINDING:
            if pos < len(seq1) and pos < len(seq2):
                res1, res2 = seq1[pos], seq2[pos]
                conserved = (res1 == res2)
                results['binding_site'].append({
                    'position': pos,
                    'DENV2': res1,
                    'DENV3': res2,
                    'conserved': conserved
                })
                if conserved:
                    results['binding_conserved_count'] += 1
        
        total_binding = len(results['binding_site'])
        results['binding_conservation_pct'] = (
            results['binding_conserved_count'] / total_binding * 100 
            if total_binding > 0 else 0
        )
        
        return results
    
    @staticmethod
    def print_alignment_summary(align_result: Dict, active_site_result: Dict):
        """打印比对摘要"""
        print("\n" + "="*70)
        print("序列比对结果")
        print("="*70)
        print(f"序列同源性: {align_result['identity_percentage']:.1f}%")
        print(f"比对长度: {align_result['length']} aa")
        
        print("\n催化三联体保守性:")
        for residue in active_site_result['catalytic_triad']:
            status = "✓" if residue['conserved'] else "✗"
            print(f"  {status} 位置 {residue['position']}: "
                  f"DENV2={residue['DENV2']}, DENV3={residue['DENV3']}")
        
        print(f"\n底物结合位点保守性: "
              f"{active_site_result['binding_conservation_pct']:.1f}% "
              f"({active_site_result['binding_conserved_count']}/{len(active_site_result['binding_site'])})")
        
        # 列出非保守残基
        non_conserved = [r for r in active_site_result['binding_site'] 
                        if not r['conserved']]
        if non_conserved:
            print("\n非保守的结合位点残基:")
            for res in non_conserved:
                print(f"  • 位置 {res['position']}: "
                      f"DENV2={res['DENV2']} → DENV3={res['DENV3']}")


# ============================================================================
# 第3步：AlphaFold预测（调用ColabFold API或本地）
# ============================================================================

class AlphaFoldPredictor:
    """使用AlphaFold进行结构预测"""
    
    @staticmethod
    def prepare_colabfold_input(sequence: str, output_file: str, job_name: str):
        """准备ColabFold输入文件"""
        with open(output_file, 'w') as f:
            f.write(f">{job_name}\n")
            # 每60个字符换行
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")
        print(f"✓ ColabFold输入已保存: {output_file}")
        print(f"\n请访问 https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb")
        print(f"上传 {output_file} 进行预测，或使用以下命令（如已安装本地版）:")
        print(f"\ncolabfold_batch {output_file} alphafold_output/ --num-models 1")
    
    @staticmethod
    def check_prediction_quality(pdb_file: str) -> Dict:
        """检查AlphaFold预测质量（pLDDT评分）"""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        plddt_scores = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.name == 'CA':  # pLDDT存储在B-factor列
                            plddt_scores.append(atom.bfactor)
        
        return {
            'mean_plddt': np.mean(plddt_scores),
            'min_plddt': np.min(plddt_scores),
            'max_plddt': np.max(plddt_scores),
            'scores': plddt_scores
        }


# ============================================================================
# 第4步：结构叠合与活性位点RMSD计算
# ============================================================================

class StructureSuperimposer:
    """结构叠合和RMSD分析"""
    
    @staticmethod
    def superimpose_structures(ref_pdb: str, mobile_pdb: str, 
                               output_aligned: str = None) -> Dict:
        """叠合两个结构并计算RMSD"""
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure('reference', ref_pdb)
        mobile_structure = parser.get_structure('mobile', mobile_pdb)
        
        ref_atoms = Selection.unfold_entities(ref_structure, 'A')
        mobile_atoms = Selection.unfold_entities(mobile_structure, 'A')
        
        # 只使用CA原子进行叠合
        ref_ca = [atom for atom in ref_atoms if atom.name == 'CA']
        mobile_ca = [atom for atom in mobile_atoms if atom.name == 'CA']
        
        # 确保原子数匹配
        min_len = min(len(ref_ca), len(mobile_ca))
        ref_ca = ref_ca[:min_len]
        mobile_ca = mobile_ca[:min_len]
        
        # 执行叠合
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_ca, mobile_ca)
        super_imposer.apply(mobile_structure.get_atoms())
        
        # 保存叠合后的结构
        if output_aligned:
            io = PDBIO()
            io.set_structure(mobile_structure)
            io.save(output_aligned)
            print(f"✓ 叠合结构已保存: {output_aligned}")
        
        return {
            'rmsd': super_imposer.rms,
            'ref_atoms': len(ref_ca),
            'mobile_atoms': len(mobile_ca)
        }
    
    @staticmethod
    def calculate_active_site_rmsd(ref_pdb: str, mobile_pdb: str, 
                                   residue_ids: List[int]) -> float:
        """计算活性位点特定残基的RMSD"""
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure('ref', ref_pdb)
        mobile_structure = parser.get_structure('mobile', mobile_pdb)
        
        ref_atoms = []
        mobile_atoms = []
        
        for res_id in residue_ids:
            try:
                ref_res = ref_structure[0]['A'][res_id]
                mobile_res = mobile_structure[0]['A'][res_id]
                
                # 使用CA原子
                ref_atoms.append(ref_res['CA'])
                mobile_atoms.append(mobile_res['CA'])
            except:
                continue
        
        if len(ref_atoms) < 3:
            return None
        
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, mobile_atoms)
        
        return super_imposer.rms


# ============================================================================
# 第5步：可视化和决策建议
# ============================================================================

class ResultVisualizer:
    """结果可视化"""
    
    @staticmethod
    def plot_plddt_scores(plddt_data: Dict, output_file: str = "plddt_plot.png"):
        """绘制pLDDT评分曲线"""
        scores = plddt_data['scores']
        
        plt.figure(figsize=(12, 4))
        plt.plot(scores, linewidth=1)
        plt.axhline(y=70, color='r', linestyle='--', label='可靠阈值 (70)')
        plt.axhline(y=90, color='g', linestyle='--', label='高置信度 (90)')
        plt.xlabel('残基位置')
        plt.ylabel('pLDDT评分')
        plt.title(f'AlphaFold预测质量 (平均pLDDT: {plddt_data["mean_plddt"]:.1f})')
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        print(f"✓ pLDDT图已保存: {output_file}")
    
    @staticmethod
    def generate_decision_report(sequence_identity: float, 
                                active_site_conservation: float,
                                catalytic_conserved: bool,
                                global_rmsd: float = None,
                                active_site_rmsd: float = None,
                                mean_plddt: float = None) -> str:
        """生成决策建议报告"""
        
        report = "\n" + "="*70 + "\n"
        report += "决策建议报告\n"
        report += "="*70 + "\n\n"
        
        # 评分系统
        score = 0
        reasons = []
        
        # 1. 序列同源性评估
        if sequence_identity >= 85:
            score += 2
            reasons.append(f"✓ 序列同源性高 ({sequence_identity:.1f}%)")
        elif sequence_identity >= 70:
            score += 1
            reasons.append(f"⚠ 序列同源性中等 ({sequence_identity:.1f}%)")
        else:
            reasons.append(f"✗ 序列同源性较低 ({sequence_identity:.1f}%)")
        
        # 2. 催化三联体保守性
        if catalytic_conserved:
            score += 3
            reasons.append("✓ 催化三联体完全保守")
        else:
            reasons.append("✗ 催化三联体存在突变（严重问题）")
        
        # 3. 活性位点保守性
        if active_site_conservation >= 85:
            score += 2
            reasons.append(f"✓ 活性位点高度保守 ({active_site_conservation:.1f}%)")
        elif active_site_conservation >= 70:
            score += 1
            reasons.append(f"⚠ 活性位点中度保守 ({active_site_conservation:.1f}%)")
        else:
            reasons.append(f"✗ 活性位点保守性低 ({active_site_conservation:.1f}%)")
        
        # 4. 结构RMSD（如果有）
        if global_rmsd is not None:
            if global_rmsd < 2.0:
                score += 2
                reasons.append(f"✓ 整体结构高度相似 (RMSD={global_rmsd:.2f}Å)")
            elif global_rmsd < 3.0:
                score += 1
                reasons.append(f"⚠ 整体结构相似度中等 (RMSD={global_rmsd:.2f}Å)")
            else:
                reasons.append(f"✗ 整体结构差异较大 (RMSD={global_rmsd:.2f}Å)")
        
        if active_site_rmsd is not None:
            if active_site_rmsd < 1.5:
                score += 2
                reasons.append(f"✓ 活性位点结构高度相似 (RMSD={active_site_rmsd:.2f}Å)")
            elif active_site_rmsd < 2.5:
                score += 1
                reasons.append(f"⚠ 活性位点结构相似度中等 (RMSD={active_site_rmsd:.2f}Å)")
        
        # 5. AlphaFold预测质量
        if mean_plddt is not None:
            if mean_plddt >= 90:
                score += 2
                reasons.append(f"✓ 预测质量优秀 (pLDDT={mean_plddt:.1f})")
            elif mean_plddt >= 70:
                score += 1
                reasons.append(f"⚠ 预测质量良好 (pLDDT={mean_plddt:.1f})")
            else:
                reasons.append(f"✗ 预测质量不佳 (pLDDT={mean_plddt:.1f})")
        
        # 打印评估详情
        report += "评估详情:\n"
        for reason in reasons:
            report += f"  {reason}\n"
        
        report += f"\n总评分: {score}/13\n\n"
        
        # 最终建议
        report += "推荐方案:\n"
        if score >= 9:
            report += "  → 使用DENV3实验结构进行对接\n"
            report += "  理由: 序列和结构高度保守，实验结构更可靠\n"
            recommendation = "DENV3"
        elif score >= 6:
            report += "  → 双轨验证：同时使用DENV2预测结构和DENV3实验结构\n"
            report += "  理由: 两者各有优势，需对接验证后选择\n"
            recommendation = "BOTH"
        else:
            report += "  → 优先使用DENV2 AlphaFold预测结构\n"
            report += "  理由: 血清型差异显著，需匹配活性数据来源\n"
            recommendation = "DENV2"
        
        if not catalytic_conserved:
            report += "\n  ⚠️ 警告: 催化三联体存在突变！建议实验验证抑制机制差异\n"
        
        report += "\n" + "="*70 + "\n"
        
        return report, recommendation


# ============================================================================
# 主流程
# ============================================================================

def main():
    """主执行流程"""
    
    print("\n🧬 DENV2 vs DENV3 结构分析流程启动")
    print("="*70)
    
    # 创建工作目录
    Path("sequences").mkdir(exist_ok=True)
    Path("structures").mkdir(exist_ok=True)
    Path("results").mkdir(exist_ok=True)
    
    # 步骤1: 获取序列
    print("\n[步骤1] 获取序列...")
    fetcher = DENVSequenceFetcher()
    sequences = fetcher.fetch_all()
    
    # 步骤2: 序列比对
    print("\n[步骤2] 序列比对分析...")
    aligner = SequenceAligner()
    align_result = aligner.align_sequences(
        sequences['DENV2'], 
        sequences['DENV3'],
        'DENV2', 'DENV3'
    )
    active_site_result = aligner.analyze_active_site(
        sequences['DENV2'], 
        sequences['DENV3']
    )
    aligner.print_alignment_summary(align_result, active_site_result)
    
    # 步骤3: 准备AlphaFold预测
    print("\n[步骤3] 准备AlphaFold预测...")
    predictor = AlphaFoldPredictor()
    predictor.prepare_colabfold_input(
        sequences['DENV2'][:185],  # NS3蛋白酶结构域前185个残基
        "sequences/DENV2_NS3_protease.fasta",
        "DENV2_NS3_protease"
    )
    
    print("\n⏸️  流程暂停：等待AlphaFold预测完成")
    print("请完成以下操作后继续:")
    print("  1. 使用上述链接或命令运行AlphaFold预测")
    print("  2. 将预测的PDB文件保存为 'structures/DENV2_predicted.pdb'")
    print("  3. 将DENV3实验结构保存为 'structures/DENV3_experimental.pdb'")
    print("  4. 运行: python this_script.py --continue")
    
    # 保存中间结果
    results = {
        'sequence_identity': align_result['identity_percentage'],
        'active_site_conservation': active_site_result['binding_conservation_pct'],
        'catalytic_conserved': active_site_result['triad_conserved']
    }
    
    with open('results/sequence_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n✓ 序列分析结果已保存到 results/sequence_analysis.json")


def continue_analysis():
    """继续分析（在获得结构后）"""
    
    print("\n🔬 继续结构比对分析...")
    
    # 加载序列分析结果
    with open('results/sequence_analysis.json', 'r') as f:
        seq_results = json.load(f)
    
    denv2_pdb = 'structures/DENV2_predicted.pdb'
    denv3_pdb = 'structures/DENV3_experimental.pdb'
    
    if not os.path.exists(denv2_pdb) or not os.path.exists(denv3_pdb):
        print("✗ 错误: 未找到结构文件，请确保已放置到structures/目录")
        return
    
    # 步骤4: 检查预测质量
    print("\n[步骤4] 检查AlphaFold预测质量...")
    predictor = AlphaFoldPredictor()
    plddt_data = predictor.check_prediction_quality(denv2_pdb)
    print(f"平均pLDDT: {plddt_data['mean_plddt']:.1f}")
    print(f"最低pLDDT: {plddt_data['min_plddt']:.1f}")
    print(f"最高pLDDT: {plddt_data['max_plddt']:.1f}")
    
    visualizer = ResultVisualizer()
    visualizer.plot_plddt_scores(plddt_data, 'results/plddt_scores.png')
    
    # 步骤5: 结构叠合
    print("\n[步骤5] 结构叠合分析...")
    superimposer = StructureSuperimposer()
    
    superimpose_result = superimposer.superimpose_structures(
        denv3_pdb, denv2_pdb, 
        'results/DENV2_aligned_to_DENV3.pdb'
    )
    print(f"整体RMSD: {superimpose_result['rmsd']:.2f} Å")
    
    # 计算活性位点RMSD
    active_site_rmsd = superimposer.calculate_active_site_rmsd(
        denv3_pdb, denv2_pdb,
        SequenceAligner.SUBSTRATE_BINDING
    )
    if active_site_rmsd:
        print(f"活性位点RMSD: {active_site_rmsd:.2f} Å")
    
    # 步骤6: 生成决策报告
    print("\n[步骤6] 生成决策建议...")
    report, recommendation = visualizer.generate_decision_report(
        sequence_identity=seq_results['sequence_identity'],
        active_site_conservation=seq_results['active_site_conservation'],
        catalytic_conserved=seq_results['catalytic_conserved'],
        global_rmsd=superimpose_result['rmsd'],
        active_site_rmsd=active_site_rmsd,
        mean_plddt=plddt_data['mean_plddt']
    )
    
    print(report)
    
    # 保存完整报告
    with open('results/decision_report.txt', 'w') as f:
        f.write(report)
    
    final_results = {
        **seq_results,
        'global_rmsd': superimpose_result['rmsd'],
        'active_site_rmsd': active_site_rmsd,
        'mean_plddt': plddt_data['mean_plddt'],
        'recommendation': recommendation
    }
    
    with open('results/final_analysis.json', 'w') as f:
        json.dump(final_results, f, indent=2)
    
    print("✓ 完整报告已保存到 results/decision_report.txt")
    print("✓ 分析数据已保存到 results/final_analysis.json")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == '--continue':
        continue_analysis()
    else:
        main()