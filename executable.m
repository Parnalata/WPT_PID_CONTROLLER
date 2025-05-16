clc; clear; close all;

Ltx = 100e-6;   
Lrx = 100e-6;   
M   = 20e-6;    
Ctx = 5;        
Crx = 7;        
Rload = 200;   

s = tf('s');
syms s_sym I_tx I_rx Vtx Vrx real

eq1 = Ltx * s_sym^2 * I_tx + M * s_sym^2 * I_rx + (1/Ctx) * I_tx == s_sym * Vtx;
eq2 = Lrx * s_sym^2 * I_rx + M * s_sym^2 * I_tx + Rload * s_sym * I_rx + (1/Crx) * I_rx == 0;

sol = solve([eq1, eq2], [I_tx, I_rx]);


Vrx_expr = Rload * sol.I_rx;
H_Vrx = simplify(Vrx_expr / Vtx);
H_Vrx_num = simplify(subs(H_Vrx));

[num, den] = numden(H_Vrx_num);
num_coeffs = sym2poly(expand(num));
den_coeffs = sym2poly(expand(den));
H = tf(num_coeffs, den_coeffs);


Kp = 0;
Ki = 1680107.5311;
Kd = 0;


tspan = [0 0.01];
x0 = [0; 0; 0; 0; 0]; 
r = @(t) 5*sin(2*pi*500*t);

[t, x] = ode45(@(t, x) wpt_dynamics_pid(t, x, Ltx, Lrx, M, Ctx, Crx, Rload, Kp, Ki, Kd, r), tspan, x0);


I_tx = x(:,1);
I_rx = x(:,2);
E_int = x(:,3);
Vrx = Rload .* I_rx;
V_Ctx = x(:, 4);
V_Crx = x(:, 5);

figure;
plot(t, Vrx, 'LineWidth', 2);
hold on;
yline(1, '--r', 'Reference');
xlabel('Time (s)');
ylabel('V_{rx} (V)');
title('Receiver Voltage V_{rx} with PID Control');
grid on;

info = stepinfo(Vrx, t, 1);
fprintf('\nStep Response Characteristics:\n');
fprintf('Settling Time: %.6f s\n', info.SettlingTime);
fprintf('Overshoot: %.2f %%\n', info.Overshoot);
fprintf('Rise Time: %.6f s\n', info.RiseTime);
fprintf('Peak Time: %.6f s\n', info.PeakTime);

C = pid(Kp, Ki, Kd);

L = C* H;
 
T = feedback(L, 1);

[numT, denT] = tfdata(T, 'v');

figure;
rlocus(T);
grid on;
title('Root Locus of Closed-Loop Transfer Function');

figure;
pzmap(T);
grid on;
title('Pole-Zero Map of Closed-Loop Transfer Function');

format long
poles = pole(H);
disp('Poles of the transfer function:');
disp(poles);