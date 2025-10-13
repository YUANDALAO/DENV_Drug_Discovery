#!/usr/bin/env python3
"""
快速分析REINVENT4运行结果 - Windows版本
适配Windows路径和CSV格式
"""

import pandas as pd
import glob
import os
from collections import Counter
import re

def analyze_results(csv_path):
    """分析单个结果CSV文件"""
    
    if not os.path.exists(csv_path):
        print(f"❌ 文件不存在: {csv_path}")
        return
    
    print(f"✅ 正在读取文件: {os.path.basename(csv_path)}\n")
    
    try:
        df = pd.read_csv(csv_path, sep=',')
    except Exception as e:
        print(f"❌ 读取CSV失败: {e}")
        return
    
    # ========================================
    # 1. 基本统计
    # ========================================
    print("=" * 80)
    print("📈 基本统计")
    print("=" * 80)
    print(f"总分子数: {len(df)}")
    print(f"平均Score: {df['Score'].mean():.4f} (范围: {df['Score'].min():.4f} - {df['Score'].max():.4f})")
    
    if 'DENV_Activity (raw)' in df.columns:
        print(f"平均pIC50: {df['DENV_Activity (raw)'].mean():.4f}")
        print(f"最高pIC50: {df['DENV_Activity (raw)'].max():.4f}")
        print(f"最低pIC50: {df['DENV_Activity (raw)'].min():.4f}")
    
    if 'SMILES_state' in df.columns:
        valid_count = (df['SMILES_state'] == 1).sum()
        print(f"有效分子: {valid_count} / {len(df)} ({valid_count/len(df)*100:.1f}%)")
    
    # ========================================
    # 2. 稳定性警报分析
    # ========================================
    print("\n" + "=" * 80)
    print("⚠️  稳定性警报分析")
    print("=" * 80)
    
    if 'matchting_patterns (Stability_Alerts)' in df.columns:
        # 统计触发警报的分子
        alert_column = df['matchting_patterns (Stability_Alerts)']
        alert_triggered = alert_column.apply(lambda x: str(x) != '[]' and pd.notna(x))
        triggered_count = alert_triggered.sum()
        
        print(f"触发警报的分子: {triggered_count} / {len(df)} ({triggered_count/len(df)*100:.1f}%)")
        
        # 统计警报类型
        all_alerts = []
        for alerts in alert_column:
            if pd.notna(alerts) and str(alerts) != '[]':
                try:
                    alert_list = eval(str(alerts))
                    if isinstance(alert_list, list):
                        all_alerts.extend(alert_list)
                except:
                    pass
        
        if all_alerts:
            alert_counter = Counter(all_alerts)
            print(f"\n最常见的警报类型 (Top 15):")
            for i, (alert, count) in enumerate(alert_counter.most_common(15), 1):
                print(f"  {i:2d}. {alert[:60]:60s} : {count:4d}次 ({count/len(df)*100:5.1f}%)")
            
            # 特别标记CF3
            cf3_alerts = sum(count for alert, count in alert_counter.items() if 'F' in alert)
            if cf3_alerts > 0:
                print(f"\n  ⚠️  含氟警报总数: {cf3_alerts}次 ({cf3_alerts/len(df)*100:.1f}%)")
        else:
            print("✅ 未发现稳定性警报")
    else:
        print("⚠️  未找到稳定性警报列")
    
    # ========================================
    # 3. R基团多样性分析
    # ========================================
    print("\n" + "=" * 80)
    print("🧬 R基团多样性分析")
    print("=" * 80)
    
    if 'R-groups' in df.columns:
        r_groups = df['R-groups'].value_counts()
        print(f"独特R基团组合数: {len(r_groups)}")
        print(f"平均每种组合出现: {len(df)/len(r_groups):.1f}次")
        
        print(f"\n前15最常见的R基团组合:")
        for i, (rg, count) in enumerate(r_groups.head(15).items(), 1):
            print(f"  {i:2d}. {str(rg)[:70]:70s} : {count:3d}次 ({count/len(df)*100:5.1f}%)")
        
        # 分析单个R基团位置
        print("\n📊 单个R基团频率分析:")
        all_r_groups = []
        for rg_str in df['R-groups']:
            if pd.notna(rg_str):
                parts = str(rg_str).split('|')
                all_r_groups.extend(parts)
        
        r_counter = Counter(all_r_groups)
        print(f"独特R基团总数: {len(r_counter)}")
        print(f"\n前20最常见的单个R基团:")
        
        cf3_total = 0
        for i, (rg, count) in enumerate(r_counter.most_common(20), 1):
            cf3_flag = "⚠️ CF3" if 'C(F)(F)F' in rg or 'C(F)F' in rg else ""
            print(f"  {i:2d}. {rg[:50]:50s} : {count:4d}次 {cf3_flag}")
            if cf3_flag:
                cf3_total += count
        
        if cf3_total > 0:
            print(f"\n  🔴 CF3基团总出现次数: {cf3_total} ({cf3_total/len(all_r_groups)*100:.1f}% of all R-groups)")
    else:
        print("⚠️  未找到R-groups列")
    
    # ========================================
    # 4. pIC50分布
    # ========================================
    print("\n" + "=" * 80)
    print("💊 pIC50活性分布")
    print("=" * 80)
    
    if 'DENV_Activity (raw)' in df.columns:
        pic50_col = df['DENV_Activity (raw)']
        
        pic50_bins = [0, 4, 5, 6, 7, 7.5, 8, 9, 15]
        pic50_labels = ['<4', '4-5', '5-6', '6-7', '7-7.5', '7.5-8', '8-9', '>9']
        df['pIC50_range'] = pd.cut(pic50_col, bins=pic50_bins, labels=pic50_labels)
        
        print("pIC50分布:")
        for label in pic50_labels:
            count = (df['pIC50_range'] == label).sum()
            bar = '█' * int(count / len(df) * 50)
            print(f"  {label:>7s}: {count:4d} ({count/len(df)*100:5.1f}%) {bar}")
        
        print(f"\n统计摘要:")
        print(f"  中位数: {pic50_col.median():.4f}")
        print(f"  标准差: {pic50_col.std():.4f}")
        print(f"  25%分位: {pic50_col.quantile(0.25):.4f}")
        print(f"  75%分位: {pic50_col.quantile(0.75):.4f}")
    else:
        print("⚠️  未找到DENV_Activity (raw)列")
    
    # ========================================
    # 5. 药物性质统计
    # ========================================
    print("\n" + "=" * 80)
    print("💊 药物性质统计 (Lipinski五规则)")
    print("=" * 80)
    
    props = {
        'MW': ('MW (raw)', 'Molecular Weight', 200, 600),
        'LogP': ('LogP (raw)', 'LogP', 0, 5),
        'HBA': ('HBA (raw)', 'H-Bond Acceptors', 0, 10),
        'HBD': ('HBD (raw)', 'H-Bond Donors', 0, 5),
        'TPSA': ('TPSA (raw)', 'TPSA', 30, 140),
        'QED': ('QED (raw)', 'Drug-likeness', 0, 1),
        'SA': ('SA (raw)', 'Synthetic Accessibility', 1, 6),
    }
    
    for name, (col, full_name, low, high) in props.items():
        if col in df.columns:
            values = df[col]
            in_range = ((values >= low) & (values <= high)).sum()
            print(f"{full_name:25s}: {values.mean():6.2f} ± {values.std():5.2f} "
                  f"(范围: {values.min():6.2f}-{values.max():6.2f}) "
                  f"[{in_range}/{len(df)} 在{low}-{high}内]")
    
    # Lipinski RO5 违规统计
    if all(col in df.columns for col in ['MW (raw)', 'LogP (raw)', 'HBA (raw)', 'HBD (raw)']):
        violations = 0
        violations += (df['MW (raw)'] > 500).sum()
        violations += (df['LogP (raw)'] > 5).sum()
        violations += (df['HBA (raw)'] > 10).sum()
        violations += (df['HBD (raw)'] > 5).sum()
        
        print(f"\n🔍 Lipinski违规: {violations}次 (平均每分子{violations/len(df):.2f}次)")
    
    # ========================================
    # 6. 导出高活性分子
    # ========================================
    print("\n" + "=" * 80)
    print("🎯 高活性分子分析")
    print("=" * 80)
    
    if 'DENV_Activity (raw)' in df.columns:
        thresholds = [7.0, 7.5, 8.0, 8.5]
        
        for threshold in thresholds:
            count = (df['DENV_Activity (raw)'] >= threshold).sum()
            print(f"pIC50 ≥ {threshold}: {count}个 ({count/len(df)*100:.1f}%)")
        
        # 取前10名
        high_activity = df.nlargest(10, 'DENV_Activity (raw)')
        
        print(f"\n🏆 Top 10 最高活性分子:")
        print("-" * 80)
        for i, (idx, row) in enumerate(high_activity.iterrows(), 1):
            print(f"\n{i}. pIC50 = {row['DENV_Activity (raw)']:.4f} | Score = {row['Score']:.4f}")
            print(f"   SMILES: {row['SMILES'][:70]}")
            if 'R-groups' in df.columns:
                print(f"   R-groups: {row['R-groups']}")
            if 'matchting_patterns (Stability_Alerts)' in df.columns:
                alerts = row['matchting_patterns (Stability_Alerts)']
                if str(alerts) != '[]':
                    print(f"   ⚠️ 警报: {alerts}")
                else:
                    print(f"   ✅ 无警报")
        
        # 保存高活性分子
        output_dir = os.path.dirname(csv_path)
        output_file = os.path.join(output_dir, "high_activity_molecules.csv")
        high_activity_all = df[df['DENV_Activity (raw)'] >= 7.5].sort_values('DENV_Activity (raw)', ascending=False)
        
        if len(high_activity_all) > 0:
            high_activity_all.to_csv(output_file, index=False)
            print(f"\n✅ {len(high_activity_all)}个高活性分子(pIC50≥7.5)已保存到:")
            print(f"   {output_file}")
    
    # ========================================
    # 7. 警告和建议
    # ========================================
    print("\n" + "=" * 80)
    print("⚡ 诊断和建议")
    print("=" * 80)
    
    warnings = []
    
    # 检查CF3
    if 'R-groups' in df.columns:
        cf3_count = df['R-groups'].astype(str).str.contains('C\(F\)\(F\)F|C\(F\)F').sum()
        cf3_ratio = cf3_count / len(df)
        if cf3_ratio > 0.5:
            warnings.append(f"🔴 CRITICAL: CF3基团过度使用 ({cf3_count}/{len(df)} = {cf3_ratio*100:.1f}%)")
            print(f"  → 建议: 添加GroupCount组件惩罚CF3")
            print(f"     weight = 1.0, smarts = 'C(F)(F)F', transform = reverse_sigmoid(high=2)")
    
    # 检查稳定性警报
    if 'matchting_patterns (Stability_Alerts)' in df.columns:
        alert_ratio = alert_triggered.sum() / len(df)
        if alert_ratio > 0.8:
            warnings.append(f"🔴 CRITICAL: 大量分子触发警报 ({alert_triggered.sum()}/{len(df)} = {alert_ratio*100:.1f}%)")
            print(f"  → 建议: 提高CustomAlerts权重从1.0到5.0")
    
    # 检查R基团多样性
    if 'R-groups' in df.columns:
        unique_rgroups = len(df['R-groups'].value_counts())
        if unique_rgroups < 50:
            warnings.append(f"🟡 WARNING: R基团多样性低 (仅{unique_rgroups}种独特组合)")
            print(f"  → 建议: 降低diversity_filter的bucket_size到4")
    
    # 检查活性
    if 'DENV_Activity (raw)' in df.columns:
        max_pic50 = df['DENV_Activity (raw)'].max()
        if max_pic50 < 8.0:
            warnings.append(f"🟡 WARNING: 未找到高活性分子 (最高pIC50={max_pic50:.2f})")
            print(f"  → 建议: 检查QSAR模型路径和权重设置")
        
        # 检查活性提升
        if len(df) > 100:
            recent_pic50 = df.tail(100)['DENV_Activity (raw)'].mean()
            early_pic50 = df.head(100)['DENV_Activity (raw)'].mean()
            improvement = recent_pic50 - early_pic50
            
            if improvement < 0.1:
                warnings.append(f"🟡 WARNING: 活性提升缓慢 (最近100 vs 最早100: +{improvement:.3f})")
                print(f"  → 建议: 检查学习率或增加QSAR权重")
    
    if not warnings:
        print("✅ 未发现严重问题")
    else:
        print(f"\n发现 {len(warnings)} 个问题需要关注")
    
    # ========================================
    # 8. 输出总结
    # ========================================
    print("\n" + "=" * 80)
    print("📋 分析总结")
    print("=" * 80)
    print(f"文件: {os.path.basename(csv_path)}")
    print(f"分子数: {len(df)}")
    if 'DENV_Activity (raw)' in df.columns:
        print(f"pIC50范围: {df['DENV_Activity (raw)'].min():.2f} - {df['DENV_Activity (raw)'].max():.2f}")
        print(f"高活性分子(≥7.5): {(df['DENV_Activity (raw)'] >= 7.5).sum()}")
    print(f"评分范围: {df['Score'].min():.4f} - {df['Score'].max():.4f}")

if __name__ == "__main__":
    import sys
    
    print("🔍 REINVENT4 结果分析器 - Windows版")
    print("=" * 80)
    
    if len(sys.argv) > 1:
        csv_path = sys.argv[1]
    else:
        # 默认路径（使用原始斜杠）
        csv_path = r"C:\Users\ucsaheu\python_projects\DENV_Drug_Discovery\05_Generative_AI_REINVENT\REINVENT4-main\experiments\runs\run9_t1200\results_1.csv"
    
    # 处理Windows路径
    csv_path = csv_path.replace('"', '').replace("'", '')
    
    analyze_results(csv_path)
    print("\n✅ 分析完成！")
    print("\n按任意键退出...")
    input()