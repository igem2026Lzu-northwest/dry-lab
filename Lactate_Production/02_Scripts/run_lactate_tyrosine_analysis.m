%% ==========================================
% 数据生成：乳酸工程 & 酪氨酸工程
% ==========================================
clear; clc; restoredefaultpath;
addpath(genpath('C:\SysBio\RAVEN'));
addpath('C:\gurobi1302\win64\matlab');
javaaddpath('C:\gurobi1302\win64\lib\gurobi.jar');

load('C:\Users\ZQQL\Desktop\igem2\01_Models\iJO1366star_melanin_final.mat');
findRxn = @(m, id) find(strcmp(m.rxns, id), 1);
ib = findRxn(model, 'BIOMASS_Ec_iJO1366_WT_53p95M');
outDir = 'C:\Users\ZQQL\Desktop\igem2\05_处理数据';
if ~exist(outDir, 'dir'), mkdir(outDir); end

% ★ 底层 Gurobi 求解器（避免 getCurrentTask 报错）
gurobi_solve = @(m) gurobi(struct('A', sparse(m.S), 'obj', full(m.c), 'lb', m.lb, 'ub', m.ub, 'rhs', m.b, 'sense', repmat('=', size(m.S,1), 1), 'modelsense', 'max', 'vtype', repmat('C', size(m.S,2), 1)), struct('OutputFlag', 0));

%% 1. 乳酸敲除测试
disp('Running Lactate Knockouts...');
m_base = model;
m_base.lb(findRxn(m_base,'EX_o2_e')) = 0; m_base.ub(findRxn(m_base,'EX_o2_e')) = 0;
m_base.c(:) = 0; m_base.c(ib) = 1;

koList = {'WT', {}; 'ΔpflB', {'PFL','OBTFL'}; 'ΔadhE', {'ACALD_TG_forward','ACALD_TG_reverse','ALCD2x_TG_forward','ALCD2x_TG_reverse'}; 'ΔackA/Δpta', {'ACKr_TG_forward','ACKr_TG_reverse','PTAr_TG_forward','PTAr_TG_reverse'}; 'Best Combo', {'PFL','OBTFL','ACALD_TG_forward','ACALD_TG_reverse','ALCD2x_TG_forward','ALCD2x_TG_reverse','FRD2','FRD3','SUCCt2_2pp'}};
exRxns = {'EX_lac__D_e','EX_for_e','EX_ac_e','EX_etoh_e'};
exNames = {'Lactate','Formate','Acetate','Ethanol'};
results = {};

for i = 1:size(koList, 1)
    m = m_base;
    for k = 1:length(koList{i,2}), j = findRxn(m, koList{i,2}{k}); if ~isempty(j), m.lb(j)=0; m.ub(j)=0; end; end
    sol = gurobi_solve(m);
    row = {koList{i,1}, sol.objval};
    for k = 1:length(exRxns)
        row{end+1} = sol.x(findRxn(m, exRxns{k}));
    end
    results(end+1, :) = row;
end
writetable(cell2table(results, 'VariableNames', [{'Condition','GrowthRate'}, exNames]), fullfile(outDir, 'table_lactate_knockouts.csv'));

%% 2. 两阶段发酵 (Stage 2 零生长)
disp('Running Two-Stage...');
m_stage2 = m_base;
m_stage2.c(:) = 0; m_stage2.c(findRxn(m_stage2,'EX_lac__D_e')) = 1;
m_stage2.lb(ib) = 0; m_stage2.ub(ib) = 0;
m_stage2.lb(findRxn(m_stage2,'EX_lac__D_e')) = 0; % 取消最小下界
sol_stage2 = gurobi_solve(m_stage2);
glc_in = -sol_stage2.x(findRxn(m_stage2,'EX_glc__D_e'));
lac_out = sol_stage2.x(findRxn(m_stage2,'EX_lac__D_e'));
twoStageTable = table({'Stage 2'}, lac_out, glc_in, lac_out/glc_in, 'VariableNames', {'Stage','Lactate','Glucose','Yield'});
writetable(twoStageTable, fullfile(outDir, 'table_two_stage.csv'));

%% 3. 酪氨酸致死性与限制测试
disp('Running Tyrosine Tests...');
m_tyr = model;
m_tyr.lb(findRxn(m_tyr,'EX_o2_e')) = -20; m_tyr.ub(findRxn(m_tyr,'EX_o2_e')) = 1000;
m_tyr.lb(findRxn(m_tyr,'ATPM')) = 0;
m_tyr.lb(findRxn(m_tyr,'EX_tyr__L_e')) = -0.15; m_tyr.ub(findRxn(m_tyr,'EX_tyr__L_e')) = 0; % 限制 Tyr
m_tyr.lb(findRxn(m_tyr,'EX_phe__L_e')) = -10; % 允许 Phe 摄取
m_tyr.lb(findRxn(m_tyr,'EX_melanin_e')) = 0.1; % 强制黑素
m_tyr.c(:) = 0; m_tyr.c(ib) = 1;

tyrConds = {'WT', {}; 'ΔpheA', {'PPND','PPNDH','HPPPNDO','PPPNDO'}; 'OE_TYRTA', {'TYRTA_TG_forward','TYRTA_TG_reverse'}; 'ΔpheA_OE_TYRTA', {'PPND','PPNDH','HPPPNDO','PPPNDO'}};
tyrResults = {};

for i = 1:size(tyrConds, 1)
    m = m_tyr;
    for k = 1:length(tyrConds{i,2}), j = findRxn(m, tyrConds{i,2}{k}); if ~isempty(j), m.lb(j)=0; m.ub(j)=0; end; end
    sol = gurobi_solve(m);
    if isfield(sol, 'x') && ~isempty(sol.x) && strcmp(sol.status, 'OPTIMAL')
        tyrResults(end+1, :) = {tyrConds{i,1}, sol.objval, sol.x(findRxn(m,'EX_tyr__L_e'))};
    else
        tyrResults(end+1, :) = {tyrConds{i,1}, 0, 'Lethal'};
    end
end
writetable(cell2table(tyrResults, 'VariableNames', {'Condition','GrowthRate','Tyr_Net_Flux'}), fullfile(outDir, 'table_tyrosine_growth.csv'));

disp('✅ 所有乳酸与酪氨酸数据 CSV 已成功生成到 05_处理数据 文件夹！');

% ★ 底层 Gurobi 求解器（避免 getCurrentTask 报错）
gurobi_solve = @(m) gurobi(struct('A', sparse(m.S), 'obj', full(m.c), 'lb', m.lb, 'ub', m.ub, 'rhs', m.b, 'sense', repmat('=', size(m.S,1), 1), 'modelsense', 'max', 'vtype', repmat('C', size(m.S,2), 1)), struct('OutputFlag', 0));

%% 1. 乳酸敲除测试
disp('Running Lactate Knockouts...');
m_base = model;
m_base.lb(findRxn(m_base,'EX_o2_e')) = 0; m_base.ub(findRxn(m_base,'EX_o2_e')) = 0;
m_base.c(:) = 0; m_base.c(ib) = 1;

koList = {'WT', {}; 'ΔpflB', {'PFL','OBTFL'}; 'ΔadhE', {'ACALD_TG_forward','ACALD_TG_reverse','ALCD2x_TG_forward','ALCD2x_TG_reverse'}; 'ΔackA/Δpta', {'ACKr_TG_forward','ACKr_TG_reverse','PTAr_TG_forward','PTAr_TG_reverse'}; 'Best Combo', {'PFL','OBTFL','ACALD_TG_forward','ACALD_TG_reverse','ALCD2x_TG_forward','ALCD2x_TG_reverse','FRD2','FRD3','SUCCt2_2pp'}};
exRxns = {'EX_lac__D_e','EX_for_e','EX_ac_e','EX_etoh_e'};
exNames = {'Lactate','Formate','Acetate','Ethanol'};
results = {};

for i = 1:size(koList, 1)
    m = m_base;
    for k = 1:length(koList{i,2}), j = findRxn(m, koList{i,2}{k}); if ~isempty(j), m.lb(j)=0; m.ub(j)=0; end; end
    sol = gurobi_solve(m);
    row = {koList{i,1}, sol.objval};
    for k = 1:length(exRxns)
        row{end+1} = sol.x(findRxn(m, exRxns{k}));
    end
    results(end+1, :) = row;
end
writetable(cell2table(results, 'VariableNames', [{'Condition','GrowthRate'}, exNames]), fullfile(outDir, 'table_lactate_knockouts.csv'));

%% 2. 两阶段发酵 (Stage 2 零生长)
disp('Running Two-Stage...');
m_stage2 = m_base;
m_stage2.c(:) = 0; m_stage2.c(findRxn(m_stage2,'EX_lac__D_e')) = 1;
m_stage2.lb(ib) = 0; m_stage2.ub(ib) = 0;
m_stage2.lb(findRxn(m_stage2,'EX_lac__D_e')) = 0; % 取消最小下界
sol_stage2 = gurobi_solve(m_stage2);
glc_in = -sol_stage2.x(findRxn(m_stage2,'EX_glc__D_e'));
lac_out = sol_stage2.x(findRxn(m_stage2,'EX_lac__D_e'));
twoStageTable = table({'Stage 2'}, lac_out, glc_in, lac_out/glc_in, 'VariableNames', {'Stage','Lactate','Glucose','Yield'});
writetable(twoStageTable, fullfile(outDir, 'table_two_stage.csv'));

%% 3. 酪氨酸致死性与限制测试
disp('Running Tyrosine Tests...');
m_tyr = model;
m_tyr.lb(findRxn(m_tyr,'EX_o2_e')) = -20; m_tyr.ub(findRxn(m_tyr,'EX_o2_e')) = 1000;
m_tyr.lb(findRxn(m_tyr,'ATPM')) = 0;
m_tyr.lb(findRxn(m_tyr,'EX_tyr__L_e')) = -0.15; m_tyr.ub(findRxn(m_tyr,'EX_tyr__L_e')) = 0; % 限制 Tyr
m_tyr.lb(findRxn(m_tyr,'EX_phe__L_e')) = -10; % 允许 Phe 摄取
m_tyr.lb(findRxn(m_tyr,'EX_melanin_e')) = 0.1; % 强制黑素
m_tyr.c(:) = 0; m_tyr.c(ib) = 1;

tyrConds = {'WT', {}; 'ΔpheA', {'PPND','PPNDH','HPPPNDO','PPPNDO'}; 'OE_TYRTA', {'TYRTA_TG_forward','TYRTA_TG_reverse'}; 'ΔpheA_OE_TYRTA', {'PPND','PPNDH','HPPPNDO','PPPNDO'}};
tyrResults = {};

for i = 1:size(tyrConds, 1)
    m = m_tyr;
    for k = 1:length(tyrConds{i,2}), j = findRxn(m, tyrConds{i,2}{k}); if ~isempty(j), m.lb(j)=0; m.ub(j)=0; end; end
    sol = gurobi_solve(m);
    if isfield(sol, 'x') && ~isempty(sol.x) && strcmp(sol.status, 'OPTIMAL')
        tyrResults(end+1, :) = {tyrConds{i,1}, sol.objval, sol.x(findRxn(m,'EX_tyr__L_e'))};
    else
        tyrResults(end+1, :) = {tyrConds{i,1}, 0, 'Lethal'};
    end
end
writetable(cell2table(results, 'VariableNames', [{'Condition','GrowthRate'}, exNames]), fullfile(outDir, 'table_lactate_knockouts.csv'));
disp('✅ 所有乳酸与酪氨酸数据 CSV 已成功生成到 05_处理数据 文件夹！');