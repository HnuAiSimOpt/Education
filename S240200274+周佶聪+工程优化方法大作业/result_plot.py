import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei']
rcParams['axes.unicode_minus'] = False

# 仅保留GA和PMA的映射
algorithm_name_map = {
    "PMA": "近端优化模因算法",
    "GA": "遗传算法"
}

# 颜色和线型配置（仅保留GA和PMA）
algorithm_style = {
    "PMA": {"color": "r", "linestyle": "-", "marker": "o"},
    "GA": {"color": "b", "linestyle": "--", "marker": "s"}
}

def load_results(algorithm_name, n_jobs=50):
    """从JSON文件加载算法结果"""
    filename = f"{algorithm_name}_results_{n_jobs}jobs.json"
    with open(filename, 'r', encoding='utf-8') as f:
        return json.load(f)

# 仅加载GA和PMA的结果
algorithm_data = {name: load_results(name) for name in algorithm_name_map.keys()}

# 1. 绘制收敛曲线对比图
plt.figure(figsize=(14, 6))

# 完工时间收敛曲线
plt.subplot(1, 2, 1)
for algo, data in algorithm_data.items():
    plt.plot(data['convergence']['generations'],
             data['convergence']['makespan'],
             color=algorithm_style[algo]["color"],
             linestyle=algorithm_style[algo]["linestyle"],
             label=algorithm_name_map[algo])
plt.xlabel('进化代数')
plt.ylabel('完工时间（分钟）')
plt.title('GA与PMA完工时间收敛曲线对比')
plt.grid(True, alpha=0.3)
plt.legend()

# 能耗收敛曲线
plt.subplot(1, 2, 2)
for algo, data in algorithm_data.items():
    plt.plot(data['convergence']['generations'],
             data['convergence']['energy'],
             color=algorithm_style[algo]["color"],
             linestyle=algorithm_style[algo]["linestyle"],
             label=algorithm_name_map[algo])
plt.xlabel('进化代数')
plt.ylabel('能耗（千瓦时）')
plt.title('GA与PMA能耗收敛曲线对比')
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig('GA_PMA_算法收敛曲线对比.png', dpi=300, bbox_inches='tight')
plt.show()

# 2. 绘制帕累托前沿对比图
plt.figure(figsize=(10, 7))

for algo, data in algorithm_data.items():
    plt.scatter(data['pareto_front']['makespan'],
                data['pareto_front']['energy'],
                c=algorithm_style[algo]["color"],
                marker=algorithm_style[algo]["marker"],
                edgecolors='k',
                alpha=0.7,
                label=algorithm_name_map[algo])

plt.xlabel('完工时间（分钟）')
plt.ylabel('能耗（千瓦时）')
plt.title('GA与PMA帕累托前沿分布对比（50个工件）')
plt.grid(True, alpha=0.3)
plt.legend()

# 标注最优解
pma_data = algorithm_data["PMA"]
min_energy = min(pma_data['pareto_front']['energy'])
min_makespan = min(pma_data['pareto_front']['makespan'])
plt.annotate('最优能耗解\n(能耗最低)',
             xy=(max(pma_data['pareto_front']['makespan']), min_energy),
             xytext=(30, -30),
             textcoords='offset points',
             arrowprops=dict(arrowstyle="->", color='k'))
plt.annotate('最快完工解\n(时间最短)',
             xy=(min_makespan, max(pma_data['pareto_front']['energy'])),
             xytext=(-80, 30),
             textcoords='offset points',
             arrowprops=dict(arrowstyle="->", color='k'))

plt.savefig('GA_PMA_帕累托前沿对比.png', dpi=300, bbox_inches='tight')
plt.show()

# 3. 绘制性能对比柱状图
metrics = {
    '算法名称': [algorithm_name_map[algo] for algo in algorithm_data.keys()],
    '完工时间': [data['performance']['best_makespan'] for data in algorithm_data.values()],
    '能耗': [data['performance']['best_energy'] for data in algorithm_data.values()],
    '运行时间': [data['performance']['runtime'] for data in algorithm_data.values()]
}

x = np.arange(len(metrics['算法名称']))
width = 0.25

fig, ax = plt.subplots(figsize=(10, 6))
rects1 = ax.bar(x - width, metrics['完工时间'], width,
                label='完工时间（分钟）', color='#1f77b4', edgecolor='navy')
rects2 = ax.bar(x, metrics['能耗'], width,
                label='能耗（千瓦时）', color='#ff7f0e', edgecolor='darkred')
rects3 = ax.bar(x + width, metrics['运行时间'], width,
                label='运行时间（秒）', color='#2ca02c', edgecolor='darkgreen')

ax.set_ylabel('性能指标')
ax.set_title('GA与PMA性能指标对比（50个工件）')
ax.set_xticks(x)
ax.set_xticklabels(metrics['算法名称'])
ax.grid(True, axis='y', alpha=0.3)
ax.legend(loc='upper right')

# 添加数据标签
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.1f}',
                    xy=(rect.get_x() + rect.get_width()/2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)

plt.tight_layout()
plt.savefig('GA_PMA_算法性能对比.png', dpi=300, bbox_inches='tight')
plt.show()