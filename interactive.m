    % 參數設定
tf = 0.28;     % m
Lf = 2;        % m
Lj = 0.5;      % m
Ls = 3;        % m
theta_fs = 50; % 度
alphas = [1, 1/2, 1/3];  % alpha 的值
theta_as = linspace(-30, 30, 10000);  % degrees
m = 1.40;
n = 3.49;
p = 1.31;
q = 4.14;

figure;
hold on;

% 繪製原始黑線
for i = 1:length(alphas)
    alpha = alphas(i);

    % 計算
    Nn_max = 3 * pi + 2 + (tf / Lf) * (alpha + (1 + alpha) / sqrt(2));
    Nt_max = 2 * alpha + 15 * (tf / Lf);
    Nm_max = (pi / 2) * (1 + (tf / Lf) ^ 2);
    c1 = sind(theta_fs + theta_as);
    c2 = cosd(theta_fs + theta_as);
    c3 = sind(theta_fs + theta_as) * ((Lj / Lf) + (Ls / Lf) * cosd(theta_fs) - (1 / 2)) - cosd(theta_fs + theta_as) * sind(theta_fs * (Ls / Lf));
    Nem = Nm_max ./ abs(c3);
    Nen = Nn_max ./ abs(c1);
    Net = Nt_max ./ abs(c2);

    % 找到每個 x 對應的 y 值的最小值
    Noninteractive = min([Nem; Nen; Net]);
    
    % noninteractive 黑線
    plot(theta_as, Noninteractive, '--', 'LineWidth', 2, 'DisplayName', ['Noninteractive St = ' num2str(i)]); 
end

% 不同 alpha 值下的解曲線
for i = 1:length(alphas)
    alpha = alphas(i);

    % 計算
    Nn_max = 3 * pi + 2 + (tf / Lf) * (alpha + (1 + alpha) / sqrt(2));
    Nt_max = 2 * alpha + 15 * (tf / Lf);
    Nm_max = (pi / 2) * (1 + (tf / Lf) ^ 2);
    c1 = sind(theta_fs + theta_as);
    c2 = cosd(theta_fs + theta_as);
    c3 = sind(theta_fs + theta_as) * ((Lj / Lf) + (Ls / Lf) * cosd(theta_fs) - (1 / 2)) - cosd(theta_fs + theta_as) * sind(theta_fs * (Ls / Lf));
    Nem = Nm_max ./ abs(c3);
    Nen = Nn_max ./ abs(c1);
    Net = Nt_max ./ abs(c2);

    % 解方程
    Interactive = zeros(size(theta_as));
    for j = 1:length(theta_as)
        % 定義方程
        equation = @(Ne) (Ne / Nen(j)) ^ q + (((Ne / Nem(j)) ^ m) + (Ne / Net(j)) ^ n) ^ (1 / p) - 1;
        
        % 使用 fsolve 解方程
        Ne_initial_guess = 1;  
        Interactive(j) = fsolve(equation, Ne_initial_guess, optimoptions('fsolve','Display','off'));
    end
    if i == 1
        save('interactive_result_st1.mat', 'Interactive');
    end
    % 將不同 alpha 值下的解曲線連線
    plot(theta_as, Interactive, 'LineWidth', 2, 'DisplayName', ['Interactive St = ' num2str(i)]); 
end

xlabel('Force Angle');
ylabel('Ne');
title('Ne vs Force Angle for Different Alphas');
ylim([0 10]); 
legend('show');  % 顯示圖例
hold off;
