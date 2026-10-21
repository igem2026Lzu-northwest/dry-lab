% 酪氨酸/黑色素 FBA 模拟
clear; clc;

load('C:\Users\ZQQL\Desktop\IGEM的WIKI\酪氨酸通路\01_模型文件\ecModel_iJO1366_melanin.mat', 'model');
bioRxn = 'BIOMASS_Ec_iJO1366_WT_53p95M';

solveModel = @(m) gurobi(...
    struct('A', sparse(m.S), 'obj', full(m.c), 'lb', m.lb, 'ub', m.ub, ...
           'rhs', m.b, 'sense', repmat('=', size(m.S,1),1), ...
           'modelsense','max', 'vtype', repmat('C', size(m.S,2),1)), ...
    struct('OutputFlag',0));

% 好氧 WT
m = model;
m = setParam(m, 'lb', 'EX_o2_e', -20);
m = setParam(m, 'lb', 'EX_glc__D_e', -10);
m = setParam(m, 'obj', bioRxn, 1);
sol = solveModel(m);
fprintf('WT aerobic: mu = %.4f\n', sol.objval);

% 酪氨酸净通量
tyrIdx = find(strcmp(m.rxns, 'EX_tyr__L_e'), 1);
if ~isempty(tyrIdx)
    fprintf('Tyr exchange: %.4f\n', sol.x(tyrIdx));
end
