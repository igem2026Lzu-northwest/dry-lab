% ============================================================
% 单菌乳酸利用测试：与酸测试相同条件（宽松营养，葡萄糖关闭，乳酸添加）
% 依赖：RAVEN 工具箱 + Gurobi
% ============================================================

clear; clc;

modelDir = 'C:\Users\ZQQL\Desktop\IGEM\community_models\SteadyCom模型\原始细菌模型';

fileNames = {
    'Acidaminococcus_fermentans_DSM_20731.xml', ...
    'Akkermansia_muciniphila_ATCC_BAA_835.xml', ...
    'Bacteroides_thetaiotaomicron_VPI_5482.xml', ...
    'Bifidobacterium_longum_E18.xml', ...
    'Eubacterium_hallii_DSM_3353.xml', ...
    'Faecalibacterium_cf_prausnitzii_KLE1255.xml', ...
    'Megasphaera_elsdenii_DSM_20460.xml', ...
    'Veillonella_parvula_Te3_DSM_2008.xml', ...
    'Lactobacillus_plantarum_JDM1.xml', ...
    'Escherichia_coli_O157_H7_str_Sakai.xml', ...
    'Roseburia_intestinalis_L1_82.xml'
};
speciesNames = {'Aci','Akk','Bth','Bif','Eub','Fae','Meg','Vei','Lac','Eco','Ros'};

lacKeywords = {'EX_lac_D', 'EX_lac_L', 'EX_lac'};

solveModel = @(model) gurobi(...
    struct('A', sparse(model.S), 'obj', full(model.c), ...
           'lb', model.lb, 'ub', model.ub, 'rhs', model.b, ...
           'sense', repmat('=', size(model.S,1), 1), ...
           'modelsense', 'max', 'vtype', repmat('C', size(model.S,2), 1)), ...
    struct('OutputFlag', 0));

mu_lactate = zeros(length(fileNames), 1);

for i = 1:length(fileNames)
    model = importModel(fullfile(modelDir, fileNames{i}));
    bioIdx = find(contains(model.rxns, 'biomass', 'IgnoreCase', true));
    bioRxn = model.rxns{bioIdx(1)};
    glcCandidates = model.rxns(contains(model.rxns, 'EX_glc') | contains(model.rxns, 'EX_Dglc'));
    glcRxn = glcCandidates{1};
    o2Candidates = model.rxns(contains(model.rxns, 'EX_o2'));
    o2Rxn = '';
    if ~isempty(o2Candidates), o2Rxn = o2Candidates{1}; end

    lacRxn = '';
    for k = 1:numel(lacKeywords)
        match = model.rxns(contains(model.rxns, lacKeywords{k}) & contains(model.rxns, 'EX_'));
        if ~isempty(match), lacRxn = match{1}; break; end
    end

    model_lac = model;
    exIdx = find(contains(model_lac.rxns, 'EX_'));
    model_lac.lb(exIdx) = -1000;
    model_lac.ub(exIdx) = 1000;
    model_lac = setParam(model_lac, 'lb', glcRxn, 0);
    model_lac = setParam(model_lac, 'ub', glcRxn, 1000);
    model_lac = setParam(model_lac, 'lb', lacRxn, -5);
    model_lac = setParam(model_lac, 'ub', lacRxn, 0);
    if ~isempty(o2Rxn)
        model_lac = setParam(model_lac, 'lb', o2Rxn, 0);
        model_lac = setParam(model_lac, 'ub', o2Rxn, 0);
    end
    model_lac = setParam(model_lac, 'obj', bioRxn, 1);

    sol = solveModel(model_lac);
    mu_lactate(i) = sol.objval;
    fprintf('%s 乳酸 μ_max = %.4f\n', speciesNames{i}, mu_lactate(i));
end

resultTable = array2table(mu_lactate, 'VariableNames', {'mu_lactate'}, 'RowNames', speciesNames);
disp(resultTable);
writetable(resultTable, fullfile(modelDir, 'lactate_mu_max_recalibrated.csv'), 'WriteRowNames', true);
