% DEA (Drag Embedment Anchor) 模擬程式

% 參數設定
Lf = 2;        % m
Ls = 3;        % m
Lj = 0.5;      % m
Af = 6;        % m^2
tf = 0.28;     % m
Wf = 3;        % m
dl = 0.073;    % m
theta_fs = 50; % 度
Su = 2;        % kPa
Su0 = 2;       % kPa
k = 0;         % kPa/m
alpha = 1;
Ncl = 12;
Nn_max = 11.98;
Nt_max = 4.39;
Nm_max = 1.645;
m = 1.56;
n = 4.19;
p = 1.57;
q = 4.43;
Za = 0;       % Zai+1 = Zai + deltaZi 記得改
theta_a = 1e-10;
theta_s = 0;
theta_as = 1e-10; % theta_as = theta_a - theta_s 記得改
En = 1;
theta_0 = 0;
theta_f = 50;
delta_t = 0.01;
last_Ne = 0;
delta_x = zeros(10000, 1);
delta_z = zeros(10000, 1);
x_positions = zeros(10000, 1);
z_positions = zeros(10000, 1);
i = 1;

% 載入互動結果
% run('interactive.m');
load('interactive_result_st1.mat','Interactive');

% 初始化第一個位置
x_positions(1) = 0;
z_positions(1) = 0;

% 主循環
while theta_f ~= 0 && i <= 10000
    fprintf('---------------------------------------------------------------------\n');
    fprintf('這是第 %d 次迴圈 \n', i);
    fprintf('theta_f = %.4f, theta_a = %.4f, Za = %.4f\n',  theta_f, theta_a, Za);
    
    % 找出相對應 desire index
    theta_as_temp = theta_as + 30;
    desire_index = round((theta_as_temp / 60) * 9999 + 1);
    desire_index = max(1, min(10000, desire_index));  % 確保索引在有效範圍內
    Ne = Interactive(desire_index);
    fprintf('Ne = %f\n', Ne);
    
    % 計算 C1, C2, C3
    C1 = sind(theta_as + theta_fs);
    C2 = cosd(theta_as + theta_fs);
    C3 = (sind(theta_fs + theta_as) * (Lj / Lf + (Ls / Lf) * cosd(theta_fs) - 0.5) - cosd(theta_fs + theta_as) * sind(theta_fs) * Ls / Lf);
    Nm = abs(Ne * C3);
    Nt = abs(Ne * C2);
    Nn = abs(Ne * C1);
    
    fprintf('C1 = %f, C2 = %f, C3 = %f\n', C1, C2, C3);
    fprintf('Nm = %f, Nt = %f, Nn = %f\n', Nm, Nt, Nn);
    
    % 計算 delta_n
    epsilon = 1e-10;
    denominator = ((Nm / Nm_max)^m + (Nt / Nt_max)^n)^((1/p) - 1);
    if denominator < epsilon || isnan(denominator) || isinf(denominator)
        delta_n = 0;
    else
        delta_n = (((Nt_max/Nn_max) * (p * q / n)) / denominator) * (((Nn / Nn_max)^(q-1)) / ((Nt / Nt_max)^(n-1)));
    end
    
    fprintf('delta_n = %f\n', delta_n);
    
    % 計算 delta_x 和 delta_z
    delta_x(i) = delta_t * cosd(theta_f) + delta_n * sind(theta_f);
    delta_z(i) = delta_t * sind(theta_f) + delta_n * cosd(theta_f);
    
    % 更新累積位置
    x_positions(i) = sum(delta_x(1:i));
    z_positions(i) = sum(delta_z(1:i));
    
    fprintf('delta_x(%d) = %f, delta_z(%d) = %f\n', i, delta_x(i), i, delta_z(i));
    
    % 計算 delta_beta
    delta_beta = ((C3 * m * Nt_max * (Nm / Nm_max)^(m-1)) / (abs(C3) * n * Nm_max * (Nt / Nt_max)^(n-1))) * (delta_t / Lf);
    fprintf('delta_beta = %f\n', delta_beta);
    
    % 計算其他參數
    delta_hat_z = delta_z(i) / dl;
    hat_Ta = Ne * Af / (dl ^ 2);
    eta = dl * k / Su;
    hat_Z = Za / dl;
    delta_Ne = abs(Ne - last_Ne);
    if i == 1
        delta_Ne = 0;
    end
    delta_hat_Ta = delta_Ne * Af / (dl ^ 2);
    delta_theta_s = delta_beta;
    
    % 檢查參數有效性
    if any(~isfinite([delta_hat_z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s]))
        warning('Some input parameters are not finite. Skipping this iteration.');
        continue;
    end
    fprintf('parameters\n');
    fprintf('delta_hat_z = %f\nEn = %f\nNcl = %f\nhat_Ta = %f\ntheta_a = %f\neta = %f\nhat_Z = %f\ndelta_hat_Ta = %f\ndelta_theta_s = %f\n', delta_hat_z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s);
    % 求解 delta_theta_a
    initial_guess = 1;
    try
        delta_theta_a = solve_delta_theta_a_complex(delta_hat_z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s, initial_guess);
    catch
        warning('Failed to solve for delta_theta_a. Skipping this iteration.');
        continue;
    end
    
    % 檢查 delta_theta_a 的有效性
    if isnan(delta_theta_a) || ~isreal(delta_theta_a) || abs(delta_theta_a) > pi
        warning('Invalid delta_theta_a calculated: %f. Capping it.', delta_theta_a);
%         delta_theta_a = sign(delta_theta_a) * min(abs(delta_theta_a), pi);
    end
    
    fprintf('delta_theta_a = %f\n', delta_theta_a);
    
    % 更新參數
    theta_a = theta_a + delta_theta_a;
    theta_f = theta_f - delta_beta;
    theta_as = theta_a - theta_s;
    theta_s = theta_s + delta_beta;
    Za = Za + delta_z(i);      
    i = i + 1;
    last_Ne = Ne;
end

% 輸出結果
if i > 10000
    disp('達到最大迭代次數');
elseif theta_f <= 0
    fprintf('theta_f = %f\n',theta_f);
    disp('計算完成：theta_f 達到或小於 0');
end

fprintf('迭代次數: %d\n', i-1);
if i > 1
    fprintf('最終 x 位置: %.4f\n', x_positions(i-1));
    fprintf('最終 z 位置: %.4f\n', z_positions(i-1));
else
    fprintf('沒有有效的計算結果\n');
end

% 繪圖
if i > 1
    figure;
    hold on;
    plot(x_positions(1:i-1), -z_positions(1:i-1), 'b-', 'LineWidth', 2);
    xlabel('拖曳距離 (Drag Distance)');
    ylabel('貫入深度 (Penetration Depth)');
    title('嵌入式拖錨 DEA 貫入軌跡示意圖');
    set(gca, 'XAxisLocation', 'top');
    grid on;
    max_x = max(x_positions(1:i-1));
    max_z = max(z_positions(1:i-1));
    xlim([0, max(max_x, 1)]);
    ylim([0, max(max_z, 1)]);
    hold off;
else
    disp('沒有足夠的數據來繪圖');
end

% 輔助函數：求解 delta_theta_a
function delta_theta_a = solve_delta_theta_a_complex(delta_hat_Z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s, initial_guess)
    equation = @(delta_theta_a) safe_equation(delta_theta_a, delta_hat_Z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s);
    options = optimset('Display', 'off');
    try
        delta_theta_a = fzero(equation, initial_guess, options);
    catch
        warning('fzero failed to find a solution. Returning initial guess.');
        delta_theta_a = initial_guess;
    end
end

function result = safe_equation(delta_theta_a, delta_hat_Z, En, Ncl, hat_Ta, theta_a, eta, hat_Z, delta_hat_Ta, delta_theta_s)
    epsilon = 1e-10;
    
    if abs(delta_hat_Z) < epsilon
        result = 0;
        return;
    end
    
    numerator = En * Ncl / hat_Ta - (((theta_a) ^ 2) / 2) * (eta + hat_Z);
    
    denominator_term = delta_theta_a - delta_theta_s;
    if abs(denominator_term) < epsilon
        result = 0;
        return;
    end
    
    denominator = theta_a + (1 / hat_Ta) * (delta_hat_Ta / denominator_term) * (((theta_a) ^ 2) / 2) * (1 - (delta_theta_s / delta_theta_a));
    
    if abs(denominator) < epsilon
        result = 0;
        return;
    end
    
    result = delta_theta_a / delta_hat_Z - numerator / denominator;
end