"""
一键运行完整诊断流程
AD分析 + 不确定性评估 + 化学空间可视化 + 综合报告
"""

import os
import sys
import pandas as pd
import json
from typing import List

# 导入前面定义的类
# from ad_analysis import ApplicabilityDomainAnalyzer
# from uncertainty_evaluation import UncertaintyEvaluator
# from chemical_space_visualization import ChemicalSpaceVisualizer
# from generate_diagnostic_report import DiagnosticReportGenerator


def run_full_diagnostic(
    training_smiles_file: str,
    generated_smiles_file: str,
    qsar_model_paths: List[str],
    output_dir: str = "diagnostic_results",
    project_name: str = "REINVENT4 Run"
):
    """
    运行完整的诊断流程
    
    Args:
        training_smiles_file: 训练集SMILES文件 (CSV, 包含'SMILES'列)
        generated_smiles_file: 生成分子文件 (CSV, 包含'SMILES'列)
        qsar_model_paths: QSAR模型文件路径列表
        output_dir: 输出目录
        project_name: 项目名称
    """
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*80)
    print(f"{project_name} - 完整诊断流程")
    print("="*80)
    print()
    
    # ========================================================================
    # Step 1: 读取数据
    # ========================================================================
    print("=" * 80)
    print("Step 1/5: 读取数据")
    print("=" * 80)
    
    try:
        train_df = pd.read_csv(training_smiles_file)
        training_smiles = train_df['SMILES'].tolist()
        print(f"✓ 训练集: {len(training_smiles)} 个分子")
        
        gen_df = pd.read_csv(generated_smiles_file)
        generated_smiles = gen_df['SMILES'].tolist()
        print(f"✓ 生成集: {len(generated_smiles)} 个分子")
        
    except Exception as e:
        print(f"✗ 读取数据失败: {e}")
        sys.exit(1)
    
    print()
    
    # ========================================================================
    # Step 2: AD分析
    # ========================================================================
    print("=" * 80)
    print("Step 2/5: 适用域(AD)分析")
    print("=" * 80)
    
    try:
        from ad_analysis import ApplicabilityDomainAnalyzer
        
        ad_analyzer = ApplicabilityDomainAnalyzer(
            training_smiles=training_smiles,
            radius=2,
            n_bits=2048
        )
        
        ad_results = ad_analyzer.calculate_ad_distance(generated_smiles)
        
        # 保存结果
        ad_path = os.path.join(output_dir, "ad_analysis_results.csv")
        ad_results.to_csv(ad_path, index=False)
        print(f"✓ AD分析完成，结果保存至: {ad_path}")
        
        # 可视化
        plot_path = os.path.join(output_dir, "ad_distribution.png")
        ad_analyzer.plot_ad_distribution(ad_results, save_path=plot_path)
        
        # 标记高风险
        high_risk = ad_analyzer.flag_high_risk_molecules(ad_results, threshold=0.4)
        high_risk_path = os.path.join(output_dir, "ad_high_risk_molecules.csv")
        high_risk.to_csv(high_risk_path, index=False)
        
    except Exception as e:
        print(f"✗ AD分析失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print()
    
    # ========================================================================
    # Step 3: 不确定性评估
    # ========================================================================
    print("=" * 80)
    print("Step 3/5: 预测不确定性评估")
    print("=" * 80)
    
    try:
        from uncertainty_evaluation import UncertaintyEvaluator
        
        # 检查模型文件
        valid_models = [p for p in qsar_model_paths if os.path.exists(p)]
        if len(valid_models) == 0:
            print("⚠️  警告: 未找到任何QSAR模型文件")
            print("   跳过不确定性评估...")
            uncertainty_results = None
        else:
            print(f"✓ 找到 {len(valid_models)} 个模型文件")
            
            evaluator = UncertaintyEvaluator(model_paths=valid_models)
            
            uncertainty_results = evaluator.predict_with_uncertainty(generated_smiles)
            
            # 保存结果
            unc_path = os.path.join(output_dir, "uncertainty_analysis.csv")
            uncertainty_results.to_csv(unc_path, index=False)
            print(f"✓ 不确定性分析完成，结果保存至: {unc_path}")
            
            # 可视化
            plot_path = os.path.join(output_dir, "uncertainty_distribution.png")
            evaluator.plot_uncertainty_analysis(uncertainty_results, save_path=plot_path)
            
            # 识别高不确定性分子
            uncertain = evaluator.identify_uncertain_predictions(uncertainty_results, threshold=1.0)
            uncertain_path = os.path.join(output_dir, "high_uncertainty_molecules.csv")
            uncertain.to_csv(uncertain_path, index=False)
            
            # 潜在遗珠
            gems = evaluator.flag_potential_gems(uncertainty_results, 
                                                 min_predicted=7.0,
                                                 min_uncertainty=0.5)
            gems_path = os.path.join(output_dir, "potential_gems.csv")
            gems.to_csv(gems_path, index=False)
        
    except Exception as e:
        print(f"✗ 不确定性评估失败: {e}")
        import traceback
        traceback.print_exc()
        uncertainty_results = None
    
    print()
    
    # ========================================================================
    # Step 4: 化学空间可视化
    # ========================================================================
    print("=" * 80)
    print("Step 4/5: 化学空间可视化")
    print("=" * 80)
    
    try:
        from chemical_space_visualization import ChemicalSpaceVisualizer
        
        # 随机采样(如果数据太大)
        max_samples = 5000
        if len(training_smiles) > max_samples:
            import random
            train_sample = random.sample(training_smiles, max_samples)
            print(f"⚠️  训练集过大，随机采样 {max_samples} 个分子")
        else:
            train_sample = training_smiles
        
        if len(generated_smiles) > max_samples:
            import random
            gen_sample = random.sample(generated_smiles, max_samples)
            print(f"⚠️  生成集过大，随机采样 {max_samples} 个分子")
        else:
            gen_sample = generated_smiles
        
        # t-SNE降维
        visualizer = ChemicalSpaceVisualizer(method='tsne')
        coords_df, reducer = visualizer.fit_transform(
            training_smiles=train_sample,
            generated_smiles=gen_sample,
            perplexity=30
        )
        
        # 保存坐标
        coords_path = os.path.join(output_dir, "chemical_space_coords.csv")
        coords_df.to_csv(coords_path, index=False)
        print(f"✓ 化学空间坐标保存至: {coords_path}")
        
        # 可视化
        plot_path = os.path.join(output_dir, "chemical_space_visualization.png")
        
        # 如果有AD和不确定性结果，合并可视化
        if uncertainty_results is not None:
            # 只保留采样的分子
            ad_sample = ad_results[ad_results['SMILES'].isin(gen_sample)]
            unc_sample = uncertainty_results[uncertainty_results['SMILES'].isin(gen_sample)]
        else:
            ad_sample = ad_results[ad_results['SMILES'].isin(gen_sample)]
            unc_sample = None
        
        visualizer.plot_chemical_space(
            df=coords_df,
            ad_results=ad_sample,
            uncertainty_results=unc_sample,
            save_path=plot_path
        )
        
        # 计算覆盖指标
        coverage_metrics = visualizer.calculate_coverage_metrics(coords_df)
        
        # 保存指标
        metrics_path = os.path.join(output_dir, "coverage_metrics.json")
        with open(metrics_path, 'w') as f:
            json.dump(coverage_metrics, f, indent=2)
        print(f"✓ 覆盖指标保存至: {metrics_path}")
        
    except Exception as e:
        print(f"✗ 化学空间可视化失败: {e}")
        import traceback
        traceback.print_exc()
        coverage_metrics = {
            'avg_nearest_neighbor_distance': 0,
            'coverage_ratio': 0,
            'extrapolation_ratio': 0,
            'n_training': len(training_smiles),
            'n_generated': len(generated_smiles)
        }
    
    print()
    
    # ========================================================================
    # Step 5: 生成综合报告
    # ========================================================================
    print("=" * 80)
    print("Step 5/5: 生成综合诊断报告")
    print("=" * 80)
    
    try:
        from generate_diagnostic_report import DiagnosticReportGenerator
        
        # 如果没有不确定性结果，创建虚拟列
        if uncertainty_results is None:
            uncertainty_results = ad_results[['SMILES']].copy()
            uncertainty_results['Mean_pIC50'] = 0.0
            uncertainty_results['Std_pIC50'] = 0.0
            uncertainty_results['Confidence'] = 'Unknown'
            print("⚠️  未进行不确定性评估，报告将只包含AD分析")
        
        reporter = DiagnosticReportGenerator(
            ad_results=ad_results,
            uncertainty_results=uncertainty_results,
            coverage_metrics=coverage_metrics,
            project_name=project_name
        )
        
        report_path = os.path.join(output_dir, "diagnostic_report.txt")
        reporter.save_report(output_path=report_path)
        
    except Exception as e:
        print(f"✗ 报告生成失败: {e}")
        import traceback
        traceback.print_exc()
    
    print()
    
    # ========================================================================
    # 完成
    # ========================================================================
    print("=" * 80)
    print("✓✓✓ 诊断完成！✓✓✓")
    print("=" * 80)
    print(f"\n所有结果保存在目录: {output_dir}/\n")
    print("主要文件:")
    print(f"  1. {output_dir}/diagnostic_report.txt - 📄 文本诊断报告")
    print(f"  2. {output_dir}/ad_distribution.png - 📊 AD分布图")
    print(f"  3. {output_dir}/uncertainty_distribution.png - 📊 不确定性分布图")
    print(f"  4. {output_dir}/chemical_space_visualization.png - 📊 化学空间图")
    print(f"  5. {output_dir}/risk_assessment.csv - 📋 风险评估表")
    print(f"  6. {output_dir}/potential_gems.csv - ⭐ 潜在遗珠列表")
    print(f"  7. {output_dir}/high_confidence_molecules.csv - ✅ 高置信度候选")
    print("=" * 80)
    print()


# ============================================================================
# 命令行接口
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="REINVENT4生成分子质量诊断工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python run_full_diagnostic.py \\
      --training data/training_set.csv \\
      --generated results/run9_t1200.csv \\
      --models model1.pkl model2.pkl model3.pkl \\
      --output diagnostic_results \\
      --name "REINVENT4 Run9_t1200"

输入文件格式:
  CSV文件必须包含'SMILES'列
  
输出:
  所有结果保存在--output指定的目录中
        """
    )
    
    parser.add_argument(
        '--training', '-t',
        required=True,
        help='训练集SMILES文件 (CSV格式，包含SMILES列)'
    )
    
    parser.add_argument(
        '--generated', '-g',
        required=True,
        help='生成分子SMILES文件 (CSV格式，包含SMILES列)'
    )
    
    parser.add_argument(
        '--models', '-m',
        nargs='+',
        default=[],
        help='QSAR模型文件路径列表 (支持多个模型用于不确定性评估)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='diagnostic_results',
        help='输出目录 (默认: diagnostic_results)'
    )
    
    parser.add_argument(
        '--name', '-n',
        default='REINVENT4 Run',
        help='项目名称 (默认: REINVENT4 Run)'
    )
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.training):
        print(f"错误: 训练集文件不存在: {args.training}")
        sys.exit(1)
    
    if not os.path.exists(args.generated):
        print(f"错误: 生成分子文件不存在: {args.generated}")
        sys.exit(1)
    
    # 运行诊断
    run_full_diagnostic(
        training_smiles_file=args.training,
        generated_smiles_file=args.generated,
        qsar_model_paths=args.models,
        output_dir=args.output,
        project_name=args.name
    )