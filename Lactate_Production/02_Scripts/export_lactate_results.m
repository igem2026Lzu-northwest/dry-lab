% 把乳酸生产 FBA 结果导出为 CSV
clear; clc;

% 假设已经有 sol 和 ecModel 在 workspace
if ~exist('sol','var') || ~exist('ecModel','var')
    error('请先运行 run_wt_anaerobic_fba 或类似脚本');
end

exch = {'EX_glc__D_e','EX_o2_e','EX_co2_e','EX_for_e','EX_succ_e', ...
        'EX_etoh_e','EX_ac_e','EX_lac__D_e','EX_h2_e'};
fluxes = zeros(numel(exch),1);
for i = 1:numel(exch)
    idx = find(strcmp(ecModel.rxns, exch{i}), 1);
    if ~isempty(idx), fluxes(i) = sol.x(idx); end
end

T = table(exch(:), fluxes, 'VariableNames', {'exchange_reaction','flux'});
writetable(T, 'lactate_flux_export.csv');
disp('已保存 lactate_flux_export.csv');
