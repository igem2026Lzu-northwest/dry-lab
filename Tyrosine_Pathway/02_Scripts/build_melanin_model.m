function community_model_improved
% ================== 11-species community model with corrected lactate parameters ==================
clear; clc;

outDir = 'C:\Users\ZQQL\Desktop\igem2';
if ~exist(outDir, 'dir'), mkdir(outDir); end

species = {'Aci','Akk','Bth','Bif','Eub','Fae','Meg','Vei','Lac','Eco','Ros'};

%% Initial conditions
X0_healthy = [0.0000, 0.0596, 0.0895, 0.0479, 0.0140, 0.0000, 0.0000, 0.0001, 0.1983, 0.0003, 0.0222];
X0_disease = [0.0001, 0.0538, 0.2142, 0.0037, 0.0207, 0.0000, 0.0001, 0.0001, 0.0488, 0.0000, 0.0191];
X0_healthy(X0_healthy == 0) = 1e-4;
X0_disease(X0_disease == 0) = 1e-4;

%% Growth parameters
mu_max_glc = [80.05, 25.0, 40.0, 80.0, 85.0, 50.13, 69.07, 107.17, 110.0, 130.0, 85.0];

% ----- 修正的乳酸 μ_max（文献值，替代宽松 FBA 结果） -----
% 顺序：Aci, Akk, Bth, Bif, Eub, Fae, Meg, Vei, Lac, Eco, Ros
mu_max_lactate_corrected = [0, 10, 5, 10, 70, 10, 70, 100, 40, 0, 50];

% 丙酸、戊酸保持之前的宽松 FBA 数据
mu_max_propionate = [62.555, 53.151, 92.502, 0, 0, 0, 62.711, 99.116, 0, 130.59, 0];
mu_max_valerate   = [79.736, 0, 104.3, 72.435, 78.93, 50.06, 69.066, 108.52, 16.249, 145.45, 69.259];

% ----- 按物种区分的 Ks（乳酸）-----
% 益生菌（Vei, Meg, Eub, Lac, Bif, Ros）Ks 低 → 优先利用
% 非益生菌（Bth, Eco, Aci, Akk, Fae）Ks 高 → 竞争劣势
Ks_lactate = [5.0, 5.0, 5.0, 0.5, 0.2, 5.0, 0.2, 0.2, 0.2, 5.0, 0.5];

Ks_glc = 0.5;
Ks_acid_default = 0.5;
Y_glc = 0.5;
Y_acid = 0.4;
D = 0.05;
S_in = 2;
sim_time = 100;
acid_prod_rate = 0.02;

% Probiotic / pathogen markers
isProbiotic = [0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1];
isPathogen  = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];

%% Healthy reference
[t_ref, Y_ref] = simulate_one_scenario(X0_healthy, 0, ...
    mu_max_glc, zeros(1,11), Ks_glc, Ks_acid_default*ones(1,11), ...
    Y_glc, Y_acid, D, S_in, sim_time);
rel_ref = Y_ref(:, 1:11) ./ sum(Y_ref(:, 1:11), 2);

%% Run three acids
acidList = {
    'Lactate',    mu_max_lactate_corrected, Ks_lactate;
    'Propionate', mu_max_propionate,        Ks_acid_default*ones(1,11);
    'Valerate',   mu_max_valerate,          Ks_acid_default*ones(1,11);
};
acidFileTag = {'lactate', 'propionate', 'valerate'};

for a = 1:size(acidList, 1)
    acidName   = acidList{a, 1};
    mu_max_a   = acidList{a, 2};
    Ks_acid_a  = acidList{a, 3};
    fileTag    = acidFileTag{a};

    fprintf('\n===== Simulating acid: %s =====\n', acidName);

    scenarios = {
        struct('name', 'Healthy',               'X0', X0_healthy, 'acid', 0);
        struct('name', 'Disease (no acid)',     'X0', X0_disease, 'acid', 0);
        struct('name', ['Disease + ' acidName], 'X0', X0_disease, 'acid', acid_prod_rate);
    };

    nScen = numel(scenarios);
    tAll = cell(nScen,1);
    relAll = cell(nScen,1);

    for s = 1:nScen
        [t, Y] = simulate_one_scenario(scenarios{s}.X0, scenarios{s}.acid, ...
            mu_max_glc, mu_max_a, Ks_glc, Ks_acid_a, ...
            Y_glc, Y_acid, D, S_in, sim_time);
        tAll{s} = t;
        relAll{s} = Y(:, 1:11) ./ sum(Y(:, 1:11), 2);
    end

    % Metrics
    probAll = cell(nScen,1);
    ppRatioAll = cell(nScen,1);
    shannonAll = cell(nScen,1);
    bcAll = cell(nScen,1);

    for s = 1:nScen
        R = relAll{s};
        probAll{s} = sum(R(:, isProbiotic==1), 2);
        proSum = sum(R(:, isProbiotic==1), 2);
        patSum = sum(R(:, isPathogen==1), 2);
        ppRatioAll{s} = proSum ./ max(patSum, 1e-6);
        shannonAll{s} = -sum(R .* log(max(R, 1e-12)), 2);
        R_ref_interp = interp1(t_ref, rel_ref, tAll{s}, 'linear', 'extrap');
        bc = zeros(numel(tAll{s}), 1);
        for i = 1:numel(tAll{s})
            bc(i) = sum(abs(R(i,:) - R_ref_interp(i,:))) / max(sum(R(i,:) + R_ref_interp(i,:)), 1e-12);
        end
        bcAll{s} = bc;
    end

    % ---- Figures ----
    colors = {'g', 'r', 'b'};
    metricNames = {'probiotic', 'ppratio', 'shannon', 'braycurtis'};
    metricData = {probAll, ppRatioAll, shannonAll, bcAll};
    metricLabels = {'Probiotic total relative abundance', ...
                    'Probiotic / Pathogen ratio', ...
                    'Shannon diversity index', ...
                    'Bray-Curtis distance to healthy'};

    for m = 1:4
        figure('Position', [100 100 900 600]); hold on;
        for s = 1:nScen
            idx = find(tAll{s} >= 1, 1);
            plot(tAll{s}(idx:end), metricData{m}{s}(idx:end), 'LineWidth', 2, 'Color', colors{s});
        end
        xlim([1 sim_time]);
        xlabel('Time (h)'); ylabel(metricLabels{m});
        legend(cellfun(@(x) x.name, scenarios, 'UniformOutput', false), 'Location', 'best');
        title(sprintf('%s with %s', metricLabels{m}, acidName));
        grid on;
        saveas(gcf, fullfile(outDir, sprintf('metric_%s_%s.png', metricNames{m}, fileTag)));
        saveas(gcf, fullfile(outDir, sprintf('metric_%s_%s.fig', metricNames{m}, fileTag)));
    end

    % 11-species dynamics
    figure('Position', [100 100 1400 900]);
    for s = 1:nScen
        subplot(1, nScen, s);
        idx = find(tAll{s} >= 1, 1);
        plot(tAll{s}(idx:end), relAll{s}(idx:end, :), 'LineWidth', 1.2);
        xlim([1 sim_time]);
        xlabel('Time (h)'); ylabel('Relative abundance');
        title(scenarios{s}.name);
        legend(species, 'Location', 'eastoutside');
        grid on;
    end
    saveas(gcf, fullfile(outDir, sprintf('all_species_%s.png', fileTag)));
    saveas(gcf, fullfile(outDir, sprintf('all_species_%s.fig', fileTag)));

    save(fullfile(outDir, sprintf('results_%s.mat', fileTag)), ...
        'tAll', 'relAll', 'probAll', 'ppRatioAll', 'shannonAll', 'bcAll', ...
        'scenarios', 'species', 'isProbiotic', 'isPathogen', 'mu_max_a', 'Ks_acid_a');
end

disp('All simulations completed.');

end

%% Helper
function [t, Y] = simulate_one_scenario(X0, acid_rate, mu_max_glc, mu_max_acid, ...
    Ks_glc, Ks_acid, Y_glc, Y_acid, D, S_in, sim_time)
    y0 = [X0, S_in, 0]';
    tspan = linspace(0, sim_time, 1000);
    [t, Y] = ode45(@(t,y) odesys(t, y, mu_max_glc, mu_max_acid, ...
        Ks_glc, Ks_acid, Y_glc, Y_acid, D, S_in, acid_rate), tspan, y0);
    [t, order] = sort(t); Y = Y(order, :);
end

%% ODE
function dydt = odesys(~, y, mu_max_glc, mu_max_acid, Ks_glc, Ks_acid, Y_glc, Y_acid, D, S_in, acid_prod_rate)
    X = y(1:11); S = y(12); A = y(13);
    mu_max_glc = mu_max_glc(:); mu_max_acid = mu_max_acid(:); Ks_acid = Ks_acid(:);
    mu = mu_max_glc .* S ./ (Ks_glc + S) + mu_max_acid .* A ./ (Ks_acid + A);
    dXdt = mu .* X - D .* X;
    dSdt = D * (S_in - S);
    for i = 1:11
        dSdt = dSdt - (mu_max_glc(i) * S / (Ks_glc + S)) * X(i) / Y_glc;
    end
    dAdt = acid_prod_rate - D * A;
    for i = 1:11
        dAdt = dAdt - (mu_max_acid(i) * A / (Ks_acid(i) + A)) * X(i) / Y_acid;
    end
    dydt = [dXdt; dSdt; dAdt];
end
