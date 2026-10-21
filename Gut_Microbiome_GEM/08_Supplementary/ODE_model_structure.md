# ODE 模型结构

## 模型类型
11 物种 Monod 动态 ODE 模型。

## 状态变量
- 各物种丰度 X_i
- 底物浓度 S_glc
- 有机酸浓度 S_acid

## 主要方程
dX_i/dt = X_i * (mu_i * S/(Ks+S) - D)
dS_glc/dt = D*(S_in - S_glc) - sum(X_i * mu_i / Y_glc)
dS_acid/dt = acid_prod_rate * sum(X_i) - D*S_acid

## 参数
见 ode_parameters.csv

## 求解
MATLAB ode45，sim_time = 100 h。
