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
theta_a = 0;
theta_s = 0;
theta_as = 0; % theta_as = theta_a - theta_s 記得改
En = 1;
theta_0 = 0;
theta_f = 50;
delta_t = 0.01;
delta_x = zeros(10000, 1);
delta_z = zeros(10000, 1);
x_positions = zeros(10000, 1);
z_positions = zeros(10000, 1);
i = 1;

% run interactive
% run('interactive.m'); %這行先註解掉避免測是浪費時間
load('interactive_result_st1.mat','Interactive');

% 初始化第一個位置
x_positions(1) = 0;
z_positions(1) = 0;

% 主循環
while theta_f ~= 0 && i <= 10000 % i 應該取10000，為了測試可調小省時間
    fprintf('---------------------------------------------------------------------\n');
    fprintf('這是第 %d 次迴圈 \n', i);
    fprintf('theta_f = %.4f, theta_a = %.4f, Za = %.4f\n',  theta_f, theta_a, Za);
    %找出相對應 desire index
    theta_as = theta_as + 30;
    desire_index = round((theta_as / 60) * 9999 + 1);
    desire_index = max(1, min(10000, desire_index));  % 確保索引在有效範圍內
    Ne = Interactive(desire_index);
    theta_as = theta_as - 30; %恢復theta_as 的值
    
    fprintf('Ne = %f\n', Ne);
    % 添加小量以避免除以零
    epsilon = 1e-10;
    fprintf('epsilon = %f\n', epsilon);
%     F = Ne * Su *Af;
%     fprintf('F = %f\n', F);
%     Fn = F * sind(theta_as + theta_fs);
%     fprintf('Fn= %f\n', Fn);
%     Ft = F * cosd(theta_as + theta_fs);
%     fprintf('Ft = %f\n', Ft);
%     M = F * Lf * (sind(theta_fs + theta_as) * (Lj / Lf + (Ls / Lf) * cosd(theta_fs) - 0.5) - cosd(theta_fs + theta_as) * sind(theta_fs) * Ls / Lf);
%     fprintf('M = %f\n', M);
%     % 檢查並處理 F 接近零的情況
%     if abs(F) < epsilon
%         C1 = 0;
%         C2 = 0;
%         C3 = 0;
%     else
%         C1 = Fn / F;
%         C2 = Ft / F;
%         C3 = M / (Lf * F);
%     end
    C1 = sind(theta_as + theta_fs);
    C2 = cosd(theta_as + theta_fs);
    C3 = (sind(theta_fs + theta_as) * (Lj / Lf + (Ls / Lf) * cosd(theta_fs) - 0.5) - cosd(theta_fs + theta_as) * sind(theta_fs) * Ls / Lf);
    Nm = abs(Ne * C3);
    Nt = abs(Ne * C2);
    Nn = abs(Ne * C1);
    fprintf('C1 = %f\n', C1);
    fprintf('C2 = %f\n', C2);
    fprintf('C3 = %f\n', C3);
    % 檢查並處理 Su, Af, 或 Lf 接近零的情況
%     if abs(Su * Af * Lf) < epsilon
%         Nm = 0;
%         Nt = 0;
%         Nn = 0;
%     else
%         Nm = abs(M / (Su * Af * Lf));
%         Nt = abs(Ft / (Su * Af));
%         Nn = abs(Fn / (Su * Af));
%     end
    fprintf('Nm = %f\n', Nm);
    fprintf('Nt = %f\n', Nt);
    fprintf('Nn = %f\n', Nn);
    % 計算 delta_n，添加檢查以避免 NaN 或 Inf
    denominator = ((Nm / Nm_max)^m + (Nt / Nt_max)^n)^((1/p) - 1);
    if denominator < epsilon || isnan(denominator) || isinf(denominator)
        delta_n = 0;
    else
        delta_n = (((Nt_max/Nn_max) * (p * q / n)) / denominator) * (((Nn / Nn_max)^(q-1)) / ((Nt / Nt_max)^(n-1)));
    end
    
    % 輸出 delta_n 的值進行調試
    fprintf('delta_n = %f\n', delta_n);
    
    delta_x(i) = delta_t * cosd(theta_f) + delta_n * sind(theta_f);
    delta_z(i) = delta_t * sind(theta_f) + delta_n * cosd(theta_f);
    
    % 更新累積位置
    x_positions(i) = sum(delta_x(1:i));
    z_positions(i) = sum(delta_z(1:i));
    
    fprintf('delta_x(%d) = %f\n', i, delta_x(i));
    fprintf('delta_z(%d) = %f\n', i, delta_z(i));
    
    % 計算 delta_beta，添加檢查以避免除以零
    if abs(C3) < epsilon || abs(Nt) < epsilon
        delta_beta = epsilon;
    else
        delta_beta = ((C3 * m * Nt_max * (Nm / Nm_max)^(m-1)) / (abs(C3) * n * Nm_max * (Nt / Nt_max)^(n-1))) * (delta_t / Lf);
    end
    fprintf('delta_beta = %f\n', delta_beta);
    
    
    hat_Ta = Ne * Af / (dl^2);
    eta = dl * k / Su;
    hat_Z = Za / dl;
    delta_hat_Z = delta_z(i) / dl;
    
    % 計算 delta_Ta，添加檢查以避免除以零
    if abs(theta_a - theta_0) < epsilon
        delta_Ta = 0;
    else
        delta_Ta = delta_z(i) * 2 * En * dl * Ncl * (Su0 + 0.5 * k * Za) / (theta_a - theta_0)^2;
    end
    fprintf('delta_Ta = %f\n', delta_Ta);
    % 檢查所有參數是否為有效數值
    if isnan(delta_z(i)) || isnan(delta_Ta) || isnan(delta_beta)
        fprintf('警告：在第 %d 次迭代中出現無效的參數值\n', i);
        fprintf('delta_z = %f, delta_Ta = %f, delta_beta = %f\n', delta_z(i), delta_Ta, delta_beta);
        break;
    end
    
    delta_theta_a = solve_equation(En, Ncl, dl, Ne, Af, theta_a, k, Su0, Za, delta_Ta, delta_beta, delta_z(i));
    
    if isnan(delta_theta_a) || delta_theta_a == 0 || abs(delta_theta_a) > pi
        fprintf('警告：在第 %d 次迭代中計算的 delta_theta_a 值無效: %f\n', i, delta_theta_a);
        break;
    end
    fprintf('delta_theta_a = %f\n', delta_theta_a);
    % 添加一個合理性檢查
    if abs(delta_theta_a) > abs(theta_f)
        fprintf('警告：在第 %d 次迭代中計算的 delta_theta_a 值過大: %f\n', i, delta_theta_a);
        delta_theta_a = sign(delta_theta_a) * min(abs(delta_theta_a), abs(theta_f));
    end

    theta_a = theta_a + delta_theta_a;
    theta_f = theta_f - delta_beta;
    theta_as = theta_a - theta_s;
    theta_s = theta_s + delta_beta;
    Za = Za + delta_z(i);      
    i = i + 1;
    
    
end

if i > 10000
    disp('達到最大迭代次數');
elseif theta_f == 0
    disp('計算完成：theta_f 達到 0');
% else
%     disp('計算因錯誤而中止');
end

% 輸出最終結果
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
    % 繪製軌跡
    plot(x_positions(1:i-1), -z_positions(1:i-1), 'b-', 'LineWidth', 2);

    % 設置軸標籤和標題
    xlabel('拖曳距離 (Drag Distance)');
    ylabel('貫入深度 (Penetration Depth)');
    title('嵌入式拖錨 DEA 貫入軌跡示意圖');

    % 移動 x 軸到頂部
    set(gca, 'XAxisLocation', 'top');
    
    grid on;

    % 設置軸的範圍（根據實際數據調整）
    max_x = max(x_positions(1:i-1));
    max_z = max(z_positions(1:i-1));
    xlim([0, max(max_x, 1)]); % 確保至少有些可見範圍
    ylim([0, max(max_z, 1)]);
    hold off;
else
    disp('失敗');
end

% helper funtion ，由於對matlab不熟，這裡主要求助於 claude
function delta_theta_a = solve_equation(En, Ncl, dl, Ne, Af, theta_a, k, Su0, Za, delta_Ta, delta_beta, delta_z)
    % 打印輸入參數
    fprintf('Input parameters:\n');
    fprintf('En=%f, Ncl=%f, dl=%f, Ne=%f, Af=%f, theta_a=%f\n', En, Ncl, dl, Ne, Af, theta_a);
    fprintf('k=%f, Su0=%f, Za=%f, delta_Ta=%f, delta_beta=%f, delta_z=%f\n', k, Su0, Za, delta_Ta, delta_beta, delta_z);
    
    % 定義方程
    fun = @(x) equation_fun(x, En, Ncl, dl, Ne, Af, theta_a, k, Su0, Za, delta_Ta, delta_beta, delta_z);
    
    % 使用較大的初始猜測值
    initial_guess = max(0.1, abs(theta_a));
    
    % 使用fsolve求解
    options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6);
    try
        [delta_theta_a, fval, exitflag, output] = fsolve(fun, initial_guess, options);
        fprintf('fsolve output: exitflag = %d, iterations = %d\n', exitflag, output.iterations);
        if exitflag <= 0
            fprintf('fsolve 無法找到解，exitflag = %d\n', exitflag);
            delta_theta_a = 0;  % 使用 0 而不是 NaN
        end
    catch e
        fprintf('Error in fsolve: %s\n', e.message);
        delta_theta_a = 0;  % 使用 0 而不是 NaN
    end
end

function F = equation_fun(x, En, Ncl, dl, Ne, Af, theta_a, k, Su0, Za, delta_Ta, delta_beta, delta_z)
    % 添加檢查以避免除以零或其他無效操作
    epsilon = 1e-10;
    
    % 左邊項
    left_term = x * dl / delta_z;
    
    % 右邊項的分子
    numerator = En * Ncl * dl / (Ne * Af) - (((theta_a) ^ 2) / 2) * (dl * k / Su0 + Za / dl);
    
    % 右邊項的分母
    denominator = theta_a;
    if abs(x - delta_beta) > epsilon && abs(Su0) > epsilon && abs(dl) > epsilon
        denominator = denominator + (dl * delta_Ta / (Ne * Af * (x - delta_beta) * Su0 * (dl) ^ 2)) * (((theta_a) ^ 2) /2) * (1 - (delta_beta / x));
    end
    
    % 計算最終的 F 值
    if abs(denominator) < epsilon
        F = left_term;
    else
        F = left_term - numerator / denominator;
    end
    
    % 檢查結果是否為有效數值
    if isnan(F) || isinf(F)
        F = realmax * sign(F);  % 使用最大的實數替代無效值
    end
end