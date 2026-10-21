% ================== 严格单菌短链脂肪酸利用测试（唯一有机碳源 = 目标酸 + CO2） ==================
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

% 目标酸：名称，交换反应ID（根据汇总表确定）
acids = {
    '甲酸',   'EX_for(e)';
    '丙酸',   'EX_ppa(e)';
    '异戊酸', 'EX_isoval(e)';   % 替代戊酸
};

% 求解函数
solveModel = @(model) gurobi(...
    struct('A', sparse(model.S), 'obj', full(model.c), ...
           'lb', model.lb, 'ub', model.ub, 'rhs', model.b, ...
           'sense', repmat('=', size(model.S,1), 1), ...
           'modelsense', 'max', 'vtype', repmat('C', size(model.S,2), 1)), ...
    struct('OutputFlag', 0));

% 结果存储：行=菌种，列=甲酸、丙酸、异戊酸
mu_results = zeros(length(fileNames), size(acids,1));
fprintf('\n===== 严格单菌短链脂肪酸利用测试 =====\n');

for sp = 1:numel(fileNames)
    model = importModel(fullfile(modelDir, fileNames{sp}));

    % 找到生物量反应
    bioIdx = find(contains(model.rxns, 'biomass', 'IgnoreCase', true));
    if isempty(bioIdx), bioIdx = find(contains(model.rxnNames, 'biomass', 'IgnoreCase', true)); end
    if isempty(bioIdx)
        warning('%s 未找到生物量反应', speciesNames{sp});
        continue;
    end
    bioRxn = model.rxns{bioIdx(1)};

    % 氧气交换
    o2Candidates = model.rxns(contains(model.rxns, 'EX_o2'));
    o2Rxn = '';
    if ~isempty(o2Candidates), o2Rxn = o2Candidates{1}; end

    % 获取所有交换反应及其代谢物化学式
    exIdx = find(contains(model.rxns, 'EX_'));
    exFormula = cell(numel(exIdx), 1);
    for j = 1:numel(exIdx)
        metIdx = find(model.S(:, exIdx(j)), 1);
        if ~isempty(metIdx) && isfield(model, 'metFormulas') && numel(model.metFormulas) >= metIdx
            exFormula{j} = model.metFormulas{metIdx};
        else
            exFormula{j} = '';
        end
    end

    % 对每种酸测试
    for a = 1:size(acids,1)
        acidName = acids{a,1};
        acidRxn = acids{a,2};

        % 检查目标酸交换反应是否存在
        if ~any(strcmp(model.rxns, acidRxn))
            warning('%s 未找到 %s 交换反应 %s，μ设为0', speciesNames{sp}, acidName, acidRxn);
            mu_results(sp,a) = 0;
            continue;
        end

        % 创建模型副本
        model_acid = model;

        % 设置每个交换反应的通量界限
        for j = 1:numel(exIdx)
            rxnIdx = exIdx(j);
            formula = exFormula{j};

            % 判断是否为无机碳源或无机营养（不含C，或为CO2）
            isInorganic = isempty(formula) || ~contains(formula, 'C') || ...
                          strcmp(formula, 'CO2') || strcmp(formula, 'CO') || ...
                          contains(formula, 'HCO3') || contains(formula, 'CO3');

            if isInorganic
                % 无机物自由交换（包括CO2）
                model_acid.lb(rxnIdx) = -1000;
                model_acid.ub(rxnIdx) = 1000;
            else
                % 有机碳源：默认关闭流入，允许流出
                model_acid.lb(rxnIdx) = 0;
                model_acid.ub(rxnIdx) = 1000;
            end
        end

        % 打开目标酸交换反应：允许流入 -5，不允许流出
        model_acid = setParam(model_acid, 'lb', acidRxn, -5);
        model_acid = setParam(model_acid, 'ub', acidRxn, 0);

        % 关闭氧气（厌氧）
        if ~isempty(o2Rxn)
            model_acid = setParam(model_acid, 'lb', o2Rxn, 0);
            model_acid = setParam(model_acid, 'ub', o2Rxn, 0);
        end

        % 设置目标函数为生物量
        model_acid = setParam(model_acid, 'obj', bioRxn, 1);

        % 求解
        sol = solveModel(model_acid);
        mu = sol.objval;
        mu_results(sp,a) = mu;

        fprintf('%s 在 %s 上的 μ_max = %.4f (反应ID: %s)\n', speciesNames{sp}, acidName, mu, acidRxn);
    end
end

% 输出表格
resultTable = array2table(mu_results, ...
    'VariableNames', {'mu_甲酸','mu_丙酸','mu_异戊酸'}, ...
    'RowNames', speciesNames);
disp(resultTable);

% 保存 CSV
outputFile = fullfile(modelDir, 'strict_acid_mu_max_results.csv');
writetable(resultTable, outputFile, 'WriteRowNames', true);
fprintf('\n结果已保存至 %s\n', outputFile);

