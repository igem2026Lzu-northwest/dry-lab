%% ==========================================
% 构建最终的 iJO1366star + Melanin 模型 (相对路径版)
% ==========================================
clear; clc;
addpath(genpath('C:\SysBio\RAVEN'));
addpath('C:\gurobi1302\win64\matlab');
javaaddpath('C:\gurobi1302\win64\lib\gurobi.jar');

% ★ 自动获取当前脚本所在目录
scriptDir = fileparts(mfilename('fullpath'));

% ★ 读取同目录下的 01_Models 里的原始 XML
xmlPath = fullfile(scriptDir, '..', '01_Models', 'iJO1366star.xml'); 
model = importModel(xmlPath);

% ★ 加上辅助函数定义！
findRxn = @(m, id) find(strcmp(m.rxns, id), 1);
findMet = @(m, id) find(strcmp(m.mets, id), 1);

% 2. 添加黑色素途径代谢物
cIdx = find(strcmp(model.comps, 'c'), 1); eIdx = find(strcmp(model.comps, 'e'), 1);
newMets = {'dopa_c', 'DOPA', 'C9H11NO4', 0, cIdx;
           'dopaquinone_c', 'DOPAquinone', 'C9H9NO4', 0, cIdx;
           'melanin_c', 'Melanin', 'C18H14N2O6', 0, cIdx;
           'melanin_e', 'Melanin', 'C18H14N2O6', 0, eIdx};
for i = 1:size(newMets, 1)
    if isempty(findMet(model, newMets{i,1}))
        model.mets{end+1} = newMets{i,1}; model.metNames{end+1} = newMets{i,2};
        model.metFormulas{end+1} = newMets{i,3}; model.metCharges(end+1) = newMets{i,4};
        model.metComps(end+1) = newMets{i,5}; model.b(end+1) = 0; model.S(end+1, :) = 0;
    end
end

% 3. 添加黑色素途径反应
rxnDefs = {
    'TYR_OX1', {'tyr__L_c', -2; 'o2_c', -1; 'dopa_c', 2}, 0, 1000;
    'TYR_OX2', {'dopa_c', -2; 'o2_c', -1; 'dopaquinone_c', 2; 'h2o_c', 2}, 0, 1000;
    'MELANIN_SYN', {'dopaquinone_c', -2; 'melanin_c', 1; 'h2o_c', 2}, 0, 1000;
    'MELANINtex', {'melanin_c', -1; 'melanin_e', 1}, 0, 1000;
    'EX_melanin_e', {'melanin_e', -1}, 0, 1000;
};
for i = 1:size(rxnDefs, 1)
    if isempty(findRxn(model, rxnDefs{i,1}))
        model.rxns{end+1} = rxnDefs{i,1}; model.rxnNames{end+1} = rxnDefs{i,1};
        model.lb(end+1) = rxnDefs{i,3}; model.ub(end+1) = rxnDefs{i,4};
        model.c(end+1) = 0; model.rev(end+1) = 0; model.grRules{end+1} = ''; model.subSystems{end+1} = {'Melanin'};
        model.S(:, end+1) = 0; colIdx = length(model.rxns);
        stoich = rxnDefs{i,2};
        for k = 1:size(stoich, 1)
            mIdx = findMet(model, stoich{k,1}); model.S(mIdx, colIdx) = stoich{k,2};
        end
    end
end

% ★ 4. 终极补丁：自动扫描并补齐所有 RAVEN 注释字段
nMets = length(model.mets);
nRxns = length(model.rxns);
fieldsToPatch = {'metMiriams', 'metNotes', 'metFrom', 'metChEBIIDs', ...
                 'rxnMiriams', 'rxnNotes', 'rxnFrom', 'eccodes', 'rxnReferences'};
for f = 1:length(fieldsToPatch)
    fieldName = fieldsToPatch{f};
    if isfield(model, fieldName)
        currLen = length(model.(fieldName));
        if currLen < nRxns && (contains(fieldName, 'rxn') || strcmp(fieldName, 'eccodes'))
            if isstruct(model.(fieldName))
                model.(fieldName)(currLen+1:nRxns) = struct('name', {}, 'value', {});
            else
                model.(fieldName)(currLen+1:nRxns) = {''};
            end
        elseif currLen < nMets && contains(fieldName, 'met')
            if isstruct(model.(fieldName))
                model.(fieldName)(currLen+1:nMets) = struct('name', {}, 'value', {});
            else
                model.(fieldName)(currLen+1:nMets) = {''};
            end
        end
    end
end

% 5. 环境与参数设置
excRxns = model.rxns(startsWith(model.rxns, 'EX_'));
for i = 1:numel(excRxns), j = findRxn(model, excRxns{i}); model.lb(j) = 0; model.ub(j) = 1000; end
whitelist = {'EX_glc__D_e','EX_o2_e','EX_nh4_e','EX_pi_e','EX_h2o_e','EX_h_e','EX_so4_e','EX_k_e','EX_fe2_e','EX_fe3_e','EX_mg2_e','EX_ca2_e','EX_mn2_e','EX_zn2_e','EX_co2_e','EX_cl_e'};
for i = 1:numel(whitelist), j = findRxn(model, whitelist{i}); if ~isempty(j), model.lb(j) = -1000; end, end
cofactors = {'EX_colipa_e','EX_thm_e','EX_thmpp_e','EX_ribflv_e','EX_fmn_e','EX_fad_e','EX_pydxn_e','EX_pydx_e','EX_nmn_e','EX_nac_e','EX_4abz_e','EX_lipoate_e','EX_gthrd_e','EX_gthox_e','EX_skm_e','EX_feenter_e','EX_enter_e','EX_mmet_e','EX_5mthf_e','EX_thf_e','EX_btn_e','EX_orot_e','EX_pnto__R_e','EX_pant__R_e','EX_adocbl_e','EX_cbl1_e'};
for i = 1:numel(cofactors), j = findRxn(model, cofactors{i}); if ~isempty(j), model.lb(j) = -0.5; end, end

model.lb(findRxn(model, 'EX_glc__D_e')) = -8;
model.lb(findRxn(model, 'ATPM')) = 8.39;
model.lb(findRxn(model, 'EX_co2_e')) = 0;

% 6. 解锁乳酸热力学拆分对
model.lb(findRxn(model, 'LDH_D_TG_forward')) = -1000; model.ub(findRxn(model, 'LDH_D_TG_forward')) = 1000;
model.lb(findRxn(model, 'LDH_D_TG_reverse')) = -1000; model.ub(findRxn(model, 'LDH_D_TG_reverse')) = 1000;

% 7. 保存与导出
outDir = fullfile(scriptDir, '..', '01_Models');
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, 'iJO1366star_melanin_final.mat'), 'model');

try
    exportModel(model, fullfile(outDir, 'iJO1366star_melanin_final.xml'));
    disp('✅ 模型已成功构建，.mat 和 .xml 均已保存！');
catch ME
    disp('⚠️ .mat 已保存，但 XML 导出仍有兼容性问题。');
end