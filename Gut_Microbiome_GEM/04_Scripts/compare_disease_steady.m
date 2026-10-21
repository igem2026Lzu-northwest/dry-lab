function plot_static_results
% 加载静态模型结果并绘制堆叠柱状图
clear; clc;

outputDir = 'C:\Users\ZQQL\Desktop\IGEM\community_models\result';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if exist('disease_compare_three_states.mat', 'file')
    load('disease_compare_three_states.mat');
else
    error('未找到 disease_compare_three_states.mat，请先运行静态模型脚本');
end

% 绘制堆叠柱状图
figure;
bar(abundanceResults', 'stacked');
set(gca, 'XTickLabel', speciesFiles(:,1));
xlabel('Species');
ylabel('Relative abundance');
legend(cellfun(@(x) x.name, scenarios, 'UniformOutput', false), 'Location', 'eastoutside');
title('Static community composition under three scenarios');
grid on;
saveas(gcf, fullfile(outputDir, 'static_three_scenarios_bar.png'));
saveas(gcf, fullfile(outputDir, 'static_three_scenarios_bar.fig'));

disp('静态模型结果图已保存。');
end