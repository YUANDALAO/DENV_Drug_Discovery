import numpy as np
import matplotlib.pyplot as plt

def sigmoid(x, low, high, k):
    mid = (high + low) / 2
    return 1 / (1 + np.exp(-k * (x - mid)))

# 当前无transform
pic50_range = np.linspace(7.0, 8.5, 100)
no_transform = pic50_range / 10

# Stage 1 transform
stage1 = sigmoid(pic50_range, low=7.0, high=8.5, k=2.0)

# Stage 2 transform
stage2 = sigmoid(pic50_range, low=7.5, high=9.0, k=2.5)

plt.figure(figsize=(12, 5))

plt.subplot(121)
plt.plot(pic50_range, no_transform, label='No Transform', linewidth=2)
plt.plot(pic50_range, stage1, label='Stage 1 (low=7.0, high=8.5, k=2.0)', linewidth=2)
plt.plot(pic50_range, stage2, label='Stage 2 (low=7.5, high=9.0, k=2.5)', linewidth=2)
plt.xlabel('pIC50', fontsize=12)
plt.ylabel('Transformed Score', fontsize=12)
plt.title('Transform Comparison', fontsize=14, fontweight='bold')
plt.legend()
plt.grid(alpha=0.3)

plt.subplot(122)
# 计算分数差异
diff_no_transform = np.diff(no_transform)
diff_stage1 = np.diff(stage1)
diff_stage2 = np.diff(stage2)

plt.plot(pic50_range[:-1], diff_no_transform, label='No Transform', linewidth=2)
plt.plot(pic50_range[:-1], diff_stage1, label='Stage 1', linewidth=2)
plt.plot(pic50_range[:-1], diff_stage2, label='Stage 2', linewidth=2)
plt.xlabel('pIC50', fontsize=12)
plt.ylabel('Score Gradient (sensitivity)', fontsize=12)
plt.title('Optimization Sensitivity', fontsize=14, fontweight='bold')
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('transform_effect.png', dpi=300)
plt.show()

# 数值对比
print("="*60)
print("Transform效果对比 (pIC50 7.5 → 8.0)")
print("="*60)
print(f"No Transform:  {7.5/10:.3f} → {8.0/10:.3f}  差距: {(8.0-7.5)/10:.3f}")
print(f"Stage 1:       {sigmoid(7.5, 7.0, 8.5, 2.0):.3f} → {sigmoid(8.0, 7.0, 8.5, 2.0):.3f}  差距: {sigmoid(8.0, 7.0, 8.5, 2.0) - sigmoid(7.5, 7.0, 8.5, 2.0):.3f}")
print(f"Stage 2:       {sigmoid(7.5, 7.5, 9.0, 2.5):.3f} → {sigmoid(8.0, 7.5, 9.0, 2.5):.3f}  差距: {sigmoid(8.0, 7.5, 9.0, 2.5) - sigmoid(7.5, 7.5, 9.0, 2.5):.3f}")