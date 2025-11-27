#!/usr/bin/env python3
"""
生成PyMOL脚本用于可视化DENV2和DENV3结构叠合
以及活性位点差异分析
"""

import json
from pathlib import Path


class PyMOLScriptGenerator:
    """生成PyMOL可视化脚本"""
    
    # 催化三联体和重要残基
    CATALYTIC_TRIAD = [51, 75, 135]
    SUBSTRATE_BINDING = [36, 37, 51, 52, 82, 84, 85, 129, 130, 132, 133, 135, 152, 155]
    
    @classmethod
    def generate_superposition_script(cls, output_file: str = "results/visualize_superposition.pml"):
        """生成结构叠合可视化脚本"""
        
        script = """# PyMOL脚本：DENV2 vs DENV3结构叠合可视化
# 使用方法：在PyMOL中运行 @visualize_superposition.pml

# 清空环境
reinitialize

# 加载结构
load structures/DENV3_experimental.pdb, DENV3
load structures/DENV2_predicted.pdb, DENV2_original
load results/DENV2_aligned_to_DENV3.pdb, DENV2_aligned

# 隐藏所有
hide everything

# 显示卡通表示
show cartoon, DENV3
show cartoon, DENV2_aligned

# 设置颜色
color cyan, DENV3
color salmon, DENV2_aligned

# 设置透明度以便比较
set cartoon_transparency, 0.3, DENV3

# 显示催化三联体（His51, Asp75, Ser135）
"""
        # 添加催化三联体显示
        for i, res_id in enumerate(cls.CATALYTIC_TRIAD):
            res_names = ['HIS', 'ASP', 'SER']
            script += f"""
# 催化残基 {res_names[i]}{res_id}
show sticks, DENV3 and resi {res_id}
show sticks, DENV2_aligned and resi {res_id}
color green, DENV3 and resi {res_id}
color yellow, DENV2_aligned and resi {res_id}
"""
        
        script += """
# 显示活性位点口袋
"""
        binding_resis = ','.join(map(str, cls.SUBSTRATE_BINDING))
        script += f"""
select active_site, resi {binding_resis}
show surface, active_site and DENV3
color lightblue, active_site and DENV3
set surface_transparency, 0.5, active_site and DENV3

# 标记关键残基
label DENV3 and resi {cls.CATALYTIC_TRIAD[0]} and name CA, "H51"
label DENV3 and resi {cls.CATALYTIC_TRIAD[1]} and name CA, "D75"
label DENV3 and resi {cls.CATALYTIC_TRIAD[2]} and name CA, "S135"

# 设置视图
zoom active_site
orient active_site

# 设置背景
bg_color white

# 创建距离标注（催化三联体之间）
distance cat_triad_1, DENV3 and resi {cls.CATALYTIC_TRIAD[0]} and name NE2, \\
                       DENV3 and resi {cls.CATALYTIC_TRIAD[1]} and name OD1
distance cat_triad_2, DENV3 and resi {cls.CATALYTIC_TRIAD[1]} and name OD2, \\
                       DENV3 and resi {cls.CATALYTIC_TRIAD[2]} and name OG

# 隐藏距离标签
hide labels, cat_triad_*

# 显示RMSD
print "\\n========== 结构叠合信息 =========="
print "绿色: DENV3实验结构"
print "粉色: DENV2预测结构（叠合后）"
print "催化三联体: His51, Asp75, Ser135"
print "================================\\n"

# 保存会话
save results/denv_superposition.pse
print "会话已保存到: results/denv_superposition.pse"
"""
        
        with open(output_file, 'w') as f:
            f.write(script)
        
        print(f"✓ PyMOL叠合脚本已保存: {output_file}")
        print(f"  使用方法: pymol {output_file}")
    
    @classmethod
    def generate_quality_check_script(cls, output_file: str = "results/check_plddt.pml"):
        """生成pLDDT质量检查脚本"""
        
        script = """# PyMOL脚本：检查AlphaFold预测质量（pLDDT着色）
# pLDDT分数存储在B-factor列中

reinitialize
load structures/DENV2_predicted.pdb, DENV2

# 按B-factor（pLDDT）着色
spectrum b, blue_white_red, DENV2, minimum=50, maximum=100

# 显示
show cartoon, DENV2
set cartoon_transparency, 0.2

# 显示关键残基
"""
        binding_resis = ','.join(map(str, cls.SUBSTRATE_BINDING))
        script += f"""
select active_site, resi {binding_resis}
show sticks, active_site
zoom active_site

# 创建pLDDT图例选区
select high_confidence, DENV2 and b > 90
select good_quality, DENV2 and b > 70 and b <= 90
select low_confidence, DENV2 and b > 50 and b <= 70
select very_low, DENV2 and b <= 50

# 打印统计信息
print "\\n========== pLDDT质量评估 =========="
print "红色区域: pLDDT > 90 (高置信度)"
print "白色区域: pLDDT 70-90 (良好)"
print "蓝色区域: pLDDT < 70 (低置信度)"
print "重点检查活性位点区域的颜色！"
print "==================================\\n"

bg_color white
save results/denv2_quality.pse
"""
        
        with open(output_file, 'w') as f:
            f.write(script)
        
        print(f"✓ PyMOL质量检查脚本已保存: {output_file}")
    
    @classmethod
    def generate_difference_analysis_script(cls, output_file: str = "results/analyze_differences.pml"):
        """生成活性位点差异分析脚本"""
        
        # 读取序列分析结果以获取非保守残基
        try:
            with open('results/sequence_analysis.json', 'r') as f:
                seq_data = json.load(f)
            
            # 这里简化处理，实际应该从更详细的数据中提取
            script = """# PyMOL脚本：活性位点差异分析

reinitialize
load structures/DENV3_experimental.pdb, DENV3
load results/DENV2_aligned_to_DENV3.pdb, DENV2

hide everything
show cartoon

# 叠合的结构
color cyan, DENV3
color salmon, DENV2
set cartoon_transparency, 0.5

# 显示所有活性位点残基
"""
            for res_id in cls.SUBSTRATE_BINDING:
                script += f"""
show sticks, resi {res_id}
"""
            
            script += """
# 如果存在非保守残基，特别标记
# （需要手动更新具体位置）
# select non_conserved, resi X+Y+Z
# color red, non_conserved
# show spheres, non_conserved and name CA

zoom resi """ + ','.join(map(str, cls.SUBSTRATE_BINDING)) + """

# 测量关键距离
# distance diff_1, DENV3 and resi X and name CA, DENV2 and resi X and name CA

bg_color white
print "检查红色标记的非保守残基！"
save results/denv_differences.pse
"""
        except:
            script = "# 请先运行主分析脚本生成序列比对数据\n"
        
        with open(output_file, 'w') as f:
            f.write(script)
        
        print(f"✓ PyMOL差异分析脚本已保存: {output_file}")
    
    @classmethod
    def generate_all_scripts(cls):
        """生成所有PyMOL脚本"""
        Path("results").mkdir(exist_ok=True)
        
        print("\n🎨 生成PyMOL可视化脚本...")
        cls.generate_superposition_script()
        cls.generate_quality_check_script()
        cls.generate_difference_analysis_script()
        
        print("\n使用说明:")
        print("1. 结构叠合查看: pymol results/visualize_superposition.pml")
        print("2. 预测质量检查: pymol results/check_plddt.pml")
        print("3. 差异分析: pymol results/analyze_differences.pml")


class HTMLReportGenerator:
    """生成交互式HTML报告"""
    
    @staticmethod
    def generate_report(output_file: str = "results/analysis_report.html"):
        """生成完整的HTML分析报告"""
        
        # 读取分析结果
        try:
            with open('results/final_analysis.json', 'r') as f:
                data = json.load(f)
        except:
            print("⚠️  未找到分析结果，请先运行主分析脚本")
            return
        
        html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DENV2 vs DENV3 结构分析报告</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .metric-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        .metric-title {{
            font-size: 14px;
            color: #666;
            margin-bottom: 5px;
        }}
        .metric-value {{
            font-size: 32px;
            font-weight: bold;
            color: #333;
        }}
        .metric-unit {{
            font-size: 18px;
            color: #999;
        }}
        .status-good {{ color: #22c55e; }}
        .status-warning {{ color: #f59e0b; }}
        .status-bad {{ color: #ef4444; }}
        .recommendation {{
            background: #fffbeb;
            border-left: 4px solid #f59e0b;
            padding: 20px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .progress-bar {{
            height: 20px;
            background: #e5e7eb;
            border-radius: 10px;
            overflow: hidden;
            margin-top: 10px;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #3b82f6, #8b5cf6);
            transition: width 0.3s ease;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e5e7eb;
        }}
        th {{
            background: #f9fafb;
            font-weight: 600;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 DENV2 vs DENV3 结构分析报告</h1>
        <p>NS2B-NS3蛋白酶比对分析与建议</p>
    </div>

    <div class="grid">
        <div class="metric-card">
            <div class="metric-title">序列同源性</div>
            <div class="metric-value {'status-good' if data.get('sequence_identity', 0) > 85 else 'status-warning' if data.get('sequence_identity', 0) > 70 else 'status-bad'}">
                {data.get('sequence_identity', 0):.1f}<span class="metric-unit">%</span>
            </div>
            <div class="progress-bar">
                <div class="progress-fill" style="width: {data.get('sequence_identity', 0)}%"></div>
            </div>
        </div>

        <div class="metric-card">
            <div class="metric-title">活性位点保守性</div>
            <div class="metric-value {'status-good' if data.get('active_site_conservation', 0) > 85 else 'status-warning' if data.get('active_site_conservation', 0) > 70 else 'status-bad'}">
                {data.get('active_site_conservation', 0):.1f}<span class="metric-unit">%</span>
            </div>
            <div class="progress-bar">
                <div class="progress-fill" style="width: {data.get('active_site_conservation', 0)}%"></div>
            </div>
        </div>

        <div class="metric-card">
            <div class="metric-title">结构RMSD</div>
            <div class="metric-value {'status-good' if data.get('global_rmsd', 999) < 2 else 'status-warning' if data.get('global_rmsd', 999) < 3 else 'status-bad'}">
                {data.get('global_rmsd', 0):.2f}<span class="metric-unit">Å</span>
            </div>
            <small>{'优秀' if data.get('global_rmsd', 999) < 2 else '良好' if data.get('global_rmsd', 999) < 3 else '偏高'}</small>
        </div>

        <div class="metric-card">
            <div class="metric-title">预测质量 (pLDDT)</div>
            <div class="metric-value {'status-good' if data.get('mean_plddt', 0) > 90 else 'status-warning' if data.get('mean_plddt', 0) > 70 else 'status-bad'}">
                {data.get('mean_plddt', 0):.1f}
            </div>
            <small>{'高置信度' if data.get('mean_plddt', 0) > 90 else '良好' if data.get('mean_plddt', 0) > 70 else '低置信'}</small>
        </div>
    </div>

    <div class="metric-card">
        <h2>📊 催化三联体检查</h2>
        <table>
            <tr>
                <th>残基位置</th>
                <th>残基类型</th>
                <th>状态</th>
            </tr>
            <tr>
                <td>His51</td>
                <td>催化组氨酸</td>
                <td class="{'status-good' if data.get('catalytic_conserved', False) else 'status-bad'}">
                    {'✓ 保守' if data.get('catalytic_conserved', False) else '✗ 突变'}
                </td>
            </tr>
            <tr>
                <td>Asp75</td>
                <td>催化天冬氨酸</td>
                <td class="{'status-good' if data.get('catalytic_conserved', False) else 'status-bad'}">
                    {'✓ 保守' if data.get('catalytic_conserved', False) else '✗ 突变'}
                </td>
            </tr>
            <tr>
                <td>Ser135</td>
                <td>催化丝氨酸</td>
                <td class="{'status-good' if data.get('catalytic_conserved', False) else 'status-bad'}">
                    {'✓ 保守' if data.get('catalytic_conserved', False) else '✗ 突变'}
                </td>
            </tr>
        </table>
    </div>

    <div class="recommendation">
        <h2>🎯 推荐方案</h2>
        <h3 style="color: #d97706; margin-top: 0;">
            {
                '方案A: 使用DENV3实验结构' if data.get('recommendation') == 'DENV3' 
                else '方案B: 双轨验证（DENV2 + DENV3）' if data.get('recommendation') == 'BOTH'
                else '方案C: 优先使用DENV2预测结构'
            }
        </h3>
        <p><strong>依据:</strong></p>
        <ul>
            {'<li>序列同源性极高，结构高度保守</li><li>实验结构质量更可靠</li>' if data.get('recommendation') == 'DENV3' else ''}
            {'<li>序列和结构存在中等差异</li><li>需要通过对接验证选择最佳方案</li><li>QSAR模型应给予较高权重</li>' if data.get('recommendation') == 'BOTH' else ''}
            {'<li>血清型差异显著</li><li>需匹配活性数据来源</li><li>主要依赖QSAR预测</li>' if data.get('recommendation') == 'DENV2' else ''}
        </ul>
        
        <h4>REINVENT4配置建议:</h4>
        <pre style="background: #1e293b; color: #e2e8f0; padding: 15px; border-radius: 5px; overflow-x: auto;">
{
    'DENV3': '''[[scoring.component]]
[scoring.component.Maize]
[[scoring.component.Maize.endpoint]]
name = "DENV3 Docking"
weight = 0.5
params.workflow = "docking_DENV3.yaml"

[[scoring.component]]
[scoring.component.QSARScorer]
[[scoring.component.QSARScorer.endpoint]]
name = "DENV2 pIC50"
weight = 0.5
params.model_path = "qsar_model.pkl"''',
    
    'BOTH': '''[[scoring.component]]
[scoring.component.QSARScorer]
[[scoring.component.QSARScorer.endpoint]]
name = "DENV2 pIC50"
weight = 0.5
params.model_path = "qsar_model.pkl"

[[scoring.component]]
[scoring.component.Maize]
[[scoring.component.Maize.endpoint]]
name = "DENV2 Docking"
weight = 0.3
params.workflow = "docking_DENV2.yaml"

[[scoring.component]]
[scoring.component.Maize]
[[scoring.component.Maize.endpoint]]
name = "DENV3 Docking"
weight = 0.2
params.workflow = "docking_DENV3.yaml"''',
    
    'DENV2': '''[[scoring.component]]
[scoring.component.QSARScorer]
[[scoring.component.QSARScorer.endpoint]]
name = "DENV2 pIC50"
weight = 0.7
params.model_path = "qsar_model.pkl"

[[scoring.component]]
[scoring.component.Maize]
[[scoring.component.Maize.endpoint]]
name = "DENV2 Docking"
weight = 0.3
params.workflow = "docking_DENV2.yaml"'''
}.get(data.get('recommendation', 'BOTH'), '')
        }</pre>
    </div>

    <div class="metric-card">
        <h2>📁 生成的文件</h2>
        <ul>
            <li><code>results/sequence_analysis.json</code> - 序列比对详细数据</li>
            <li><code>results/plddt_scores.png</code> - AlphaFold质量评估图</li>
            <li><code>results/DENV2_aligned_to_DENV3.pdb</code> - 叠合后的DENV2结构</li>
            <li><code>results/decision_report.txt</code> - 文本格式完整报告</li>
            <li><code>results/visualize_superposition.pml</code> - PyMOL可视化脚本</li>
        </ul>
    </div>

    <div class="metric-card">
        <h2>🔬 下一步操作</h2>
        <ol>
            <li><strong>结构可视化:</strong> 使用PyMOL检查叠合质量
                <pre style="background: #f3f4f6; padding: 10px; margin-top: 5px;">pymol results/visualize_superposition.pml</pre>
            </li>
            <li><strong>准备对接:</strong> 根据推荐方案配置MAIZE工作流</li>
            <li><strong>验证方案:</strong> 用5-10个已知抑制剂测试对接</li>
            <li><strong>集成REINVENT4:</strong> 使用上述TOML配置</li>
        </ol>
    </div>

    <footer style="text-align: center; color: #999; margin-top: 50px; padding: 20px;">
        <p>生成时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>DENV Structure Analysis Pipeline v1.0</p>
    </footer>

    <script>
        // 添加简单的交互效果
        document.querySelectorAll('.metric-card').forEach(card => {{
            card.addEventListener('mouseenter', function() {{
                this.style.transform = 'translateY(-5px)';
                this.style.transition = 'transform 0.3s ease';
            }});
            card.addEventListener('mouseleave', function() {{
                this.style.transform = 'translateY(0)';
            }});
        }});

        // 动画进度条
        window.addEventListener('load', function() {{
            document.querySelectorAll('.progress-fill').forEach(bar => {{
                const width = bar.style.width;
                bar.style.width = '0';
                setTimeout(() => {{
                    bar.style.width = width;
                }}, 100);
            }});
        }});
    </script>
</body>
</html>"""
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✓ HTML报告已生成: {output_file}")
        print(f"  在浏览器中打开查看: file://{Path(output_file).absolute()}")


def main():
    """生成所有可视化脚本和报告"""
    print("\n" + "="*70)
    print("DENV 结构分析 - 可视化脚本生成器")
    print("="*70)
    
    # 生成PyMOL脚本
    PyMOLScriptGenerator.generate_all_scripts()
    
    # 生成HTML报告
    print("\n📄 生成HTML交互式报告...")
    HTMLReportGenerator.generate_report()
    
    print("\n" + "="*70)
    print("✨ 所有可视化文件已生成完成!")
    print("="*70)
    print("\n推荐查看顺序:")
    print("1. 在浏览器打开: results/analysis_report.html (整体概览)")
    print("2. PyMOL查看结构: pymol results/visualize_superposition.pml")
    print("3. 检查预测质量: pymol results/check_plddt.pml")
    print("4. 详细文本报告: cat results/decision_report.txt")


if __name__ == "__main__":
    main()