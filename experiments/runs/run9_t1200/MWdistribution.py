import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("full_evaluation/all_molecules_metrics.csv")

# 按QED分组
qed_bins = [0, 0.3, 0.5, 0.7, 1.0]
qed_labels = ['Poor(<0.3)', 'Fair(0.3-0.5)', 'Good(0.5-0.7)', 'Excellent(>0.7)']

df['QED_Category'] = pd.cut(df['QED'], bins=qed_bins, labels=qed_labels)

# 箱线图
fig, ax = plt.subplots(figsize=(10, 6))
df.boxplot(column='MW', by='QED_Category', ax=ax)
ax.set_xlabel('QED Category', fontweight='bold', fontsize=12)
ax.set_ylabel('Molecular Weight', fontweight='bold', fontsize=12)
ax.set_title('MW Distribution by QED Category', fontweight='bold', fontsize=14)
plt.suptitle('')  # 移除默认标题
plt.tight_layout()
plt.savefig('mw_by_qed.png', dpi=300)
print("✓ 图表已保存: mw_by_qed.png")