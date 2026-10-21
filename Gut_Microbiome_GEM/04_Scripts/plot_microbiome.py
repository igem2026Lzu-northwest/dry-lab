import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch
import matplotlib as mpl

# ================= 1. 全局配置与 Times New Roman 字体 =================
mpl.rcParams['font.family'] = 'serif'
# 优先使用 Times New Roman，如果系统没安装则回退到 DejaVu Serif
mpl.rcParams['font.serif'] = ['Times New Roman', 'SimSun', 'DejaVu Serif']
mpl.rcParams['mathtext.fontset'] = 'stix'  # 数学公式贴近 Times 风格
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['figure.dpi'] = 300

# 自动获取当前脚本所在目录 (即 04_Scripts)
current_dir = os.path.dirname(os.path.abspath(__file__))
# 向上退一级 (Gut_Microbiome_GEM)，然后拼接 07_Figures
out_base = os.path.join(current_dir, '..', '07_Figures')
os.makedirs(out_base, exist_ok=True)

# ================= 2. 核心参数定义 =================
species = ['Aci','Akk','Bth','Bif','Eub','Fae','Meg','Vei','Lac','Eco','Ros']
# 功能群划分：Eco和Aci为机会致病菌，Bth不计入（条件共生）
isProbiotic = np.array([0,1,0,1,1,1,1,1,1,0,1], dtype=bool)
isPathogen  = np.array([1,0,0,0,0,0,0,0,0,1,0], dtype=bool)

# 16S 原始数据
X0_healthy_raw = np.array([0.0000, 0.0596, 0.0895, 0.0479, 0.0140, 0.0000, 0.0000, 0.0001, 0.1983, 0.0003, 0.0222])
X0_disease_raw = np.array([0.0001, 0.0538, 0.2142, 0.0037, 0.0207, 0.0000, 0.0001, 0.0001, 0.0488, 0.0000, 0.0191])

# ★ 修复：严格归一化 (Sum=1)
X0_healthy = X0_healthy_raw / np.sum(X0_healthy_raw)
X0_disease = X0_disease_raw / np.sum(X0_disease_raw)
# 防止 0 值导致 ODE 求解失败
X0_healthy[X0_healthy < 1e-5] = 1e-4
X0_disease[X0_disease < 1e-5] = 1e-4
X0_healthy = X0_healthy / np.sum(X0_healthy)
X0_disease = X0_disease / np.sum(X0_disease)

mu_max_glc = np.array([80.05, 25.0, 40.0, 80.0, 85.0, 50.13, 69.07, 107.17, 110.0, 130.0, 85.0])

# ★ 修复：使用真实的、未篡改的乳酸参数（与 Wiki 表格完全一致）
mu_max_acid_dict = {
    'Lactate':    np.array([78.587, 54.948, 104.3, 72.559, 78.93, 50.132, 69.066, 108.72, 68.459, 145.78, 69.373]),
    'Propionate': np.array([62.555, 53.151, 92.502, 0, 0, 0, 62.711, 99.116, 0, 130.59, 0]),
    'Valerate':   np.array([79.736, 0, 104.3, 72.435, 78.93, 50.06, 69.066, 108.52, 16.249, 145.45, 69.259])
}

Ks_glc, Ks_acid = 0.5, 0.5
Y_glc, Y_acid = 0.5, 0.4
D, S_in, sim_time, acid_prod_rate = 0.05, 2.0, 100.0, 0.02

# ================= 3. ODE 求解与模拟 =================
def odesys(y, t, mu_max_glc, mu_max_acid, acid_rate):
    X = y[:11]; S = y[11]; A = y[12]
    mu_glc = mu_max_glc * S / (Ks_glc + S)
    mu_acid = mu_max_acid * A / (Ks_acid + A)
    mu = mu_glc + mu_acid
    dXdt = mu * X - D * X
    dSdt = D * (S_in - S) - np.sum(mu_glc * X / Y_glc)
    dAdt = acid_rate - D * A - np.sum(mu_acid * X / Y_acid)
    return np.concatenate([dXdt, [dSdt, dAdt]])

def simulate(X0, acid_rate, mu_max_a, y0_steady_override=None):
    # 先跑 50h 让系统达到稳态，避免初始瞬态震荡影响比较
    if y0_steady_override is not None:
        y0_steady = y0_steady_override
    else:
        y0 = np.concatenate([X0, [S_in, 0.0]])
        t_burn = np.linspace(0, 50, 1000)
        Y_burn = odeint(odesys, y0, t_burn, args=(mu_max_glc, mu_max_a, acid_rate))
        y0_steady = np.clip(Y_burn[-1, :], 0, None)
    
    t = np.linspace(0, sim_time, 1000)
    Y = odeint(odesys, y0_steady, t, args=(mu_max_glc, mu_max_a, acid_rate))
    Y = np.clip(Y, 0, None)
    
    rel = Y[:, :11] / Y[:, :11].sum(axis=1, keepdims=True)
    return t, rel, y0_steady

# ================= 4. 配色 =================
species_colors = ['#D55E00', '#009E73', '#999999', '#56B4E9', '#66CC88', '#4C8CCC', '#A673D9', '#E69933', '#F2BF1A', '#CC0000', '#808080']

# ================= 5. 主循环生成单酸图表 =================
comparison_data = {}  # 用于最终横向对比收集数据
healthy_baseline = {'Probiotic': None, 'Ratio': None, 'Shannon': None, 'BC': None}

for acid_name, mu_max_a in mu_max_acid_dict.items():
    acid_dir = os.path.join(out_base, acid_name)
    os.makedirs(acid_dir, exist_ok=True)
    
    # 模拟三种场景
    t, rel_H, _ = simulate(X0_healthy, 0.0, mu_max_a)
    _, rel_D, steady_disease = simulate(X0_disease, 0.0, mu_max_a)
    # 加酸组：必须共享疾病稳态起点，确保对比公平
    _, rel_DA, _ = simulate(None, acid_prod_rate, mu_max_a, y0_steady_override=steady_disease)
    
    rel_dict = {'Healthy': rel_H, 'Disease': rel_D, f'Disease + {acid_name}': rel_DA}
    metrics = {}
    
    for scen, rel in rel_dict.items():
        prob = rel[:, isProbiotic].sum(axis=1)
        pat = rel[:, isPathogen].sum(axis=1)
        ratio = prob / np.maximum(pat, 1e-3)
        shannon = -np.sum(rel * np.log(np.maximum(rel, 1e-12)), axis=1)
        bc = np.abs(rel - rel_H).sum(axis=1) / np.maximum((rel + rel_H).sum(axis=1), 1e-12)
        metrics[scen] = {'Probiotic': prob, 'Ratio': ratio, 'Shannon': shannon, 'BC': bc}
        
        # 保存 CSV
        df = pd.DataFrame(rel, columns=species)
        df.insert(0, 'Time', t)
        df.to_csv(os.path.join(acid_dir, f'{scen}_dynamics.csv'), index=False)
    
    # 记录数据供汇总图使用
    comparison_data[acid_name] = {
        'Probiotic': metrics[f'Disease + {acid_name}']['Probiotic'],
        'Ratio': metrics[f'Disease + {acid_name}']['Ratio'],
        'Shannon': metrics[f'Disease + {acid_name}']['Shannon'],
        'BC': metrics[f'Disease + {acid_name}']['BC']
    }
    if healthy_baseline['Probiotic'] is None:
        for key in healthy_baseline.keys():
            healthy_baseline[key] = metrics['Healthy'][key]

    # ---------- 图 A：四指标拼图 (单酸) ----------
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    colors = ['#009E73', '#D55E00', '#0072B2']
    configs = [('Probiotic', 'Probiotic Abundance', 'Relative Abundance', axes[0,0]),
               ('Ratio', 'Probiotic / Pathogen Ratio', 'Ratio', axes[0,1]),
               ('Shannon', 'Shannon Diversity Index', 'Index', axes[1,0]),
               ('BC', 'Bray-Curtis Distance to Healthy', 'Distance', axes[1,1])]
    
    for i, (key, title, ylab, ax) in enumerate(configs):
        for j, scen in enumerate(rel_dict.keys()):
            ax.plot(t, metrics[scen][key], linewidth=2.5, color=colors[j], label=scen)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Time (h)', fontsize=12)
        ax.set_ylabel(ylab, fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(0, sim_time)
        ax.margins(x=0)
        if i == 0: ax.legend(frameon=False, fontsize=11)
    
    fig.suptitle(f'Community Metrics under {acid_name}', fontsize=16, fontweight='bold', y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(os.path.join(acid_dir, f'Metrics_{acid_name}.png'), transparent=True)
    plt.close()

    # ---------- 图 B：11物种堆叠面积图 ----------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    for idx, scen in enumerate(rel_dict.keys()):
        ax = axes[idx]
        ax.stackplot(t, rel_dict[scen].T, colors=species_colors, linewidth=0.2, edgecolor='none')
        ax.set_xlim(0, sim_time)
        ax.set_ylim(0, 1)
        ax.set_title(scen, fontsize=14, fontweight='bold')
        ax.set_xlabel('Time (h)', fontsize=12)
        ax.margins(x=0)
        if idx == 0: ax.set_ylabel('Relative Abundance', fontsize=12)
            
    legend_patches = [Patch(color=species_colors[i], label=species[i]) for i in range(11)]
    # 图例防遮挡：放到图表右侧边缘之外
    fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1.02, 0.5), 
               frameon=False, fontsize=12, title="Species", title_fontsize=14)
    fig.suptitle(f'Species Dynamics under {acid_name}', fontsize=18, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 0.90, 0.92]) 
    plt.savefig(os.path.join(acid_dir, f'Dynamics_{acid_name}.png'), bbox_inches='tight', transparent=True)
    plt.close()

print("🎉 单酸群落图表已生成完毕！")

# ================= 6. 终极对比图：健康基线 vs 三种酸干预 =================
fig, axes = plt.subplots(2, 2, figsize=(16, 10))
comp_colors = {'Lactate': '#0072B2', 'Propionate': '#E69F00', 'Valerate': '#CC79A7'}

metric_configs = [
    ('Probiotic', 'Probiotic Abundance', 'Relative Abundance', axes[0,0]),
    ('Ratio', 'Probiotic / Pathogen Ratio', 'Ratio', axes[0,1]),
    ('Shannon', 'Shannon Diversity Index', 'Index', axes[1,0]),
    ('BC', 'Bray-Curtis Distance to Healthy', 'Distance', axes[1,1])
]

for i, (key, title, ylab, ax) in enumerate(metric_configs):
    # 画健康基线（灰色粗虚线）
    ax.plot(t, healthy_baseline[key], linewidth=3.5, color='#999999', linestyle='--', label='Healthy (Baseline)')
    
    # 画三种酸干预组
    for acid_name in ['Lactate', 'Propionate', 'Valerate']:
        ax.plot(t, comparison_data[acid_name][key], linewidth=2.5, color=comp_colors[acid_name], 
                label=f'Disease + {acid_name}')
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Time (h)', fontsize=12)
    ax.set_ylabel(ylab, fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(0, sim_time)
    ax.margins(x=0)
    if i == 0:
        ax.legend(frameon=False, fontsize=11, loc='lower right')

fig.suptitle('Comprehensive Metrics: Healthy Baseline vs Acid Interventions', fontsize=18, fontweight='bold', y=0.96)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(os.path.join(out_base, 'All_Metrics_Comparison_vs_Healthy.png'), transparent=True)
plt.close()

print("🎉 终极对比图 All_Metrics_Comparison_vs_Healthy.png 已生成完毕！")