import scipy.stats as stats

# 相关性分析
correlation = df_all[['Score', 'QED', 'MW', 'LogP']].corr()

print("相关性矩阵:")
print(correlation)

# 可视化
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Score vs QED
axes[0].hexbin(df_all['Score'], df_all['QED'], gridsize=50, cmap='YlOrRd')
axes[0].set_xlabel('Score', fontweight='bold')
axes[0].set_ylabel('QED', fontweight='bold')
axes[0].set_title(f'Score vs QED (r={correlation.loc["Score", "QED"]:.3f})', 
                  fontweight='bold')

# Score vs MW
axes[1].hexbin(df_all['Score'], df_all['MW'], gridsize=50, cmap='YlGnBu')
axes[1].set_xlabel('Score', fontweight='bold')
axes[1].set_ylabel('MW', fontweight='bold')
axes[1].set_title(f'Score vs MW (r={correlation.loc["Score", "MW"]:.3f})', 
                  fontweight='bold')

# QED vs MW
axes[2].hexbin(df_all['QED'], df_all['MW'], gridsize=50, cmap='viridis')
axes[2].set_xlabel('QED', fontweight='bold')
axes[2].set_ylabel('MW', fontweight='bold')
axes[2].set_title(f'QED vs MW (r={correlation.loc["QED", "MW"]:.3f})', 
                  fontweight='bold')

plt.tight_layout()
plt.savefig('correlation_analysis.png', dpi=300)
print("✓ 相关性图表已保存: correlation_analysis.png")