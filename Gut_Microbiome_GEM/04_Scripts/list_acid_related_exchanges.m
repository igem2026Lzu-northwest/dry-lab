% ================== 列出所有菌种中酸相关及无机碳交换反应 ==================
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

acidKeywords = {'for', 'form', 'ppa', 'prop', 'val', 'isoval', 'pentanoate', 'co2', 'hco3', 'bicarbonate'};

allRelated = {};

for sp = 1:numel(fileNames)
    model = importModel(fullfile(modelDir, fileNames{sp}));
    exIdx = find(contains(model.rxns, 'EX_'));

    fprintf('\n===== 物种：%s =====\n', speciesNames{sp});
    fprintf('总交换反应数：%d\n', numel(exIdx));

    countC = 0;
    for j = 1:numel(exIdx)
        rxnIdx = exIdx(j);
        metIdx = find(model.S(:, rxnIdx), 1);
        if isempty(metIdx), continue; end
        formula = '';
        if isfield(model, 'metFormulas') && numel(model.metFormulas) >= metIdx
            formula = model.metFormulas{metIdx};
        end
        if contains(formula, 'C'), countC = countC + 1; end
    end
    fprintf('含碳交换反应数：%d，不含碳：%d\n', countC, numel(exIdx)-countC);

    related = [];
    for j = 1:numel(exIdx)
        rxnIdx = exIdx(j);
        rxnName = model.rxns{rxnIdx};
        if any(cellfun(@(k) contains(rxnName, k, 'IgnoreCase', true), acidKeywords))
            related(end+1) = rxnIdx;
        end
    end

    fprintf('酸相关/无机碳相关交换反应：\n');
    if isempty(related)
        fprintf('  无\n');
    else
        for r = 1:numel(related)
            rxnIdx = related(r);
            rxnName = model.rxns{rxnIdx};
            metIdx = find(model.S(:, rxnIdx), 1);
            metID = model.mets{metIdx};
            metName = '';
            if isfield(model, 'metNames') && numel(model.metNames) >= metIdx
                metName = model.metNames{metIdx};
            end
            formula = '';
            if isfield(model, 'metFormulas') && numel(model.metFormulas) >= metIdx
                formula = model.metFormulas{metIdx};
            end
            fprintf('  %s | %s | %s | %s\n', rxnName, metID, metName, formula);
            allRelated{end+1,1} = speciesNames{sp};
            allRelated{end,2} = rxnName;
            allRelated{end,3} = metID;
            allRelated{end,4} = formula;
        end
    end
end

if ~isempty(allRelated)
    resultTable = cell2table(allRelated, 'VariableNames', {'Species','RxnID','MetID','Formula'});
    writetable(resultTable, fullfile(modelDir, 'related_exchange_reactions_summary.csv'));
    fprintf('\n汇总表已保存至 related_exchange_reactions_summary.csv\n');
end
