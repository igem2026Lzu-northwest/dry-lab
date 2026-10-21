% ================== 构建8菌联合模型并保存 ==================
clear; clc;

global CBT_LP_SOLVER CBT_MILP_SOLVER CBT_QP_SOLVER CBT_NLP_SOLVER
CBT_LP_SOLVER = 'gurobi';
CBT_MILP_SOLVER = 'gurobi';
CBT_QP_SOLVER = 'gurobi';
CBT_NLP_SOLVER = 'gurobi';

modelDir = 'C:\Users\ZQQL\Desktop\IGEM\community_models';

speciesFiles = {
    'Aci', 'Acidaminococcus_fermentans_DSM_20731.xml', 'biomass400';
    'Akk', 'Akkermansia_muciniphila_ATCC_BAA_835.xml', 'biomass345';
    'Bth', 'Bacteroides_thetaiotaomicron_VPI_5482.xml', 'biomass345';
    'Bif', 'Bifidobacterium_longum_E18.xml', 'biomass345';
    'Eub', 'Eubacterium_hallii_DSM_3353.xml', 'biomass536';
    'Fae', 'Faecalibacterium_cf_prausnitzii_KLE1255.xml', 'biomass525';
    'Meg', 'Megasphaera_elsdenii_DSM_20460.xml', 'biomass370';
    'Vei', 'Veillonella_parvula_Te3_DSM_2008.xml', 'biomass205';
};

% 1. 加载并单独合并每个模型
jointModels = cell(8,1);
for i = 1:8
    tag = speciesFiles{i,1};
    fname = speciesFiles{i,2};
    bio = speciesFiles{i,3};
    model = readCbModel(fullfile(modelDir, fname));
    jointModels{i} = createMultipleSpeciesModel({model}, {bio}, 'nameTagsModels', {tag});
    jointModels{i} = ensureFields(jointModels{i});
end

% 2. 合并所有单独合并的模型
modelJoint = jointModels{1};
for i = 2:8
    modelJoint = mergeModels({modelJoint, jointModels{i}});
end

% 保存模型
save('modelJoint_8spp.mat', 'modelJoint');
disp('模型已保存为 modelJoint_8spp.mat');

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
