% ================== 双菌模型：乳酸梯度曲线 ==================
clear; clc;

% 设置求解器
global CBT_LP_SOLVER CBT_MILP_SOLVER CBT_QP_SOLVER CBT_NLP_SOLVER
CBT_LP_SOLVER = 'gurobi';
CBT_MILP_SOLVER = 'gurobi';
CBT_QP_SOLVER = 'gurobi';
CBT_NLP_SOLVER = 'gurobi';

modelDir = 'C:\Users\ZQQL\Desktop\IGEM\community_models';

% 读取模型
Bth = readCbModel(fullfile(modelDir, 'Bacteroides_thetaiotaomicron_VPI_5482.xml'));
Vei = readCbModel(fullfile(modelDir, 'Veillonella_parvula_Te3_DSM_2008.xml'));

% 分别合并单个模型
jointBth = createMultipleSpeciesModel({Bth}, {'biomass345'}, 'nameTagsModels', {'Bth'});
jointVei = createMultipleSpeciesModel({Vei}, {'biomass205'}, 'nameTagsModels', {'Vei'});

% 补全字段
jointBth = ensureFields(jointBth);
jointVei = ensureFields(jointVei);

% 合并两个联合模型
modelJoint = mergeModels({jointBth, jointVei});

% 生物量反应索引
idxBthBiomass = find(strcmp(modelJoint.rxns, 'Bthbiomass345'));
idxVeiBiomass = find(strcmp(modelJoint.rxns, 'Veibiomass205'));
biomassIdx = [idxBthBiomass; idxVeiBiomass];

% 扫描参数
glcLevel = 10;                % 固定葡萄糖
lacLevels = [0 1 5 10 20 50]; % 乳酸梯度

abundanceBth = zeros(size(lacLevels));
abundanceVei = zeros(size(lacLevels));

for li = 1:numel(lacLevels)
    model = modelJoint;  % 直接赋值，不是 copy()
    
    % 设置葡萄糖
    glcRxns = model.rxns(contains(model.rxns, 'EX_glc_D'));
    for i = 1:numel(glcRxns)
        model = changeRxnBounds(model, glcRxns{i}, -glcLevel, 'l');
    end
    
    % 设置乳酸
    lacRxns = model.rxns(contains(model.rxns, 'EX_lac_D') | contains(model.rxns, 'EX_lac_L'));
    for i = 1:numel(lacRxns)
        model = changeRxnBounds(model, lacRxns{i}, -lacLevels(li), 'l');
    end
    
    % 氧气为0
    o2Rxns = model.rxns(contains(model.rxns, 'EX_o2'));
    for i = 1:numel(o2Rxns)
        model = changeRxnBounds(model, o2Rxns{i}, 0, 'b');
    end
    
    % 设置每个物种最小生长率（可适当提高，如0.1，观察影响）
    model.lb(biomassIdx) = 0.01;
    model.ub(biomassIdx) = 1000;
    
    % 目标函数：最大化总生物量
    model.c(:) = 0;
    model.c(biomassIdx) = 1;
    
    % 求解
    sol = optimizeCbModel(model, 'max');
    
    if sol.stat == 1
        if isfield(sol, 'full')
            flux = sol.full;
        elseif isfield(sol, 'x')
            flux = sol.x;
        elseif isfield(sol, 'v')
            flux = sol.v;
        else
            error('找不到解向量字段');
        end
        bthFlux = flux(idxBthBiomass);
        veiFlux = flux(idxVeiBiomass);
        total = bthFlux + veiFlux;
        abundanceBth(li) = bthFlux / total;
        abundanceVei(li) = veiFlux / total;
    else
        abundanceBth(li) = NaN;
        abundanceVei(li) = NaN;
    end
end

% 绘制曲线
figure;
plot(lacLevels, abundanceBth, '-o', 'LineWidth', 2); hold on;
plot(lacLevels, abundanceVei, '-s', 'LineWidth', 2);
xlabel('Lactate uptake rate (mmol/gDW/h)');
ylabel('Relative abundance');
legend('Bth', 'Vei');
title(sprintf('Glucose uptake = %d, effect of lactate on relative abundance', glcLevel));
grid on;

% 保存数据
save('abundance_curve.mat', 'lacLevels', 'abundanceBth', 'abundanceVei');

% ================== 局部函数：ensureFields ==================
function model = ensureFields(model)
    if ~isfield(model, 'id'), model.id = 'model'; end
    if ~isfield(model, 'rev'), model.rev = double(model.lb < 0); end
    if ~isfield(model, 'metComps'), model.metComps = ones(numel(model.mets),1); end
    if ~isfield(model, 'rxnGeneMat'), model.rxnGeneMat = sparse(numel(model.rxns), numel(model.genes)); end
    if ~isfield(model, 'grRules'), model.grRules = repmat({''}, numel(model.rxns),1); end
    if ~isfield(model, 'rxnFrom'), model.rxnFrom = repmat({model.id}, numel(model.rxns),1); end
    if ~isfield(model, 'metFrom'), model.metFrom = repmat({model.id}, numel(model.mets),1); end
    if ~isfield(model, 'geneFrom'), model.geneFrom = repmat({model.id}, numel(model.genes),1); end
    if ~isfield(model, 'subSystems'), model.subSystems = repmat({{''}}, numel(model.rxns),1); end
    if ~isfield(model, 'rxnNames'), model.rxnNames = model.rxns; end
    if ~isfield(model, 'compNames'), model.compNames = model.comps; end
    if ~isfield(model, 'rxnNotes'), model.rxnNotes = repmat({''}, numel(model.rxns),1); end
    if ~isfield(model, 'metNames'), model.metNames = model.mets; end
    if ~isfield(model, 'metFormulas'), model.metFormulas = repmat({''}, numel(model.mets),1); end
    if ~isfield(model, 'metCharges'), model.metCharges = zeros(numel(model.mets),1); end
    if ~isfield(model, 'description'), model.description = ''; end
    if ~isfield(model, 'modelVersion'), model.modelVersion = '1.0'; end
end