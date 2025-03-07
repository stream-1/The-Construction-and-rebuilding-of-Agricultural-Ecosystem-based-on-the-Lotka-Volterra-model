% figure 10 and figure 11
clear,clc;
% parameter setting
r1 = 1.3; r2 = 1.6; r3 = 1.2; r4 = 0.08; r5 = 0.1; r6 = 0.4; 
K1 = 10; K2 = 25; K3 = 25; K4 = 0.25; K5 = 0.2; K6 = 10; % Amended environmental capacity from: K1 = 50; K2 = 25; K3 = 25; K4 = 3; K5 = 2;
% theta1 = 0.2;theta2 = 1;theta3 = 1;theta4 =  0.0833 ;theta5 = 0.1;
gama1 = 1.4; gama2 = 0.1;  gama3 = 0.11;  gama4 = 0.09; 
alpha1 = 0.06; alpha2 = 0.06; alpha3 = 0.1; 
beta1 = 0.8; beta2 = 0.6; beta3 = 0.5;
sigema1=0.08; sigema2=0.08; sigema3=0.08; sigema4 = 0.22; sigema5 = 0.43;

% initial conditions
C = 25; 
W = 11; 
I = 16; 
B = 2; 
R = 1; 
F = 4; 
state = [2, 2, 2, 2, 2, 2];

% time parameter
h = 0.1; 
num_steps = 365; 

% Storing Historical Values of Variables
states_history = zeros(num_steps, 6);
C_history = zeros(num_steps, 1);
W_history = zeros(num_steps, 1);
I_history = zeros(num_steps, 1);
B_history = zeros(num_steps, 1);
R_history = zeros(num_steps, 1);
F_history = zeros(num_steps, 1);

% Define the transfer matrix 
P_C = [0.7, 0.2, 0.1; 0.4, 0.5, 0.1; 0.2, 0.4, 0.4]; 
P_W = [0.6, 0.3, 0.1; 0.1, 0.5, 0.4; 0.1, 0.1, 0.8]; 
P_I = [0.8, 0.15, 0.05; 0.2, 0.6, 0.2; 0.1, 0.25, 0.65]; 
P_B = [0.8, 0.15, 0.05; 0.25, 0.5, 0.25; 0.15, 0.3, 0.55]; 
P_R = [0.8, 0.15, 0.05; 0.3, 0.5, 0.2; 0.1, 0.3, 0.6]; 
P_F = [0.75, 0.2, 0.05; 0.2, 0.6, 0.2; 0.15, 0.25, 0.6]; 

% Store all transfer matrices in the cell array
transition_matrices = {P_C, P_W, P_I, P_B, P_R, P_F};

% iterative computation
for t = 1:num_steps
    % Record status history
    states_history(t, :) = state;
    C_history(t) = C;
    W_history(t) = W;
    I_history(t) = I;
    B_history(t) = B;
    R_history(t) = R;
    F_history(t) = F;
    
    % Updating variables using the Lungkuta method
    % 1. dC/dt
k1_C = r1 * C * (1 - C / K1 - gama1 * W / K2) - alpha1 * C * I + beta1 * C + sigema4 * B * C + sigema5 * F * C;
k2_C = r1 * (C + 0.5 * h * k1_C) * (1 - (C + 0.5 * h * k1_C) / K1 - gama1 * W / K2) - ...
       alpha1 * (C + 0.5 * h * k1_C) * I + beta1 * (C + 0.5 * h * k1_C) + sigema4 * B * (C + 0.5 * h * k1_C) + ...
       sigema5 * F * (C + 0.5 * h * k1_C);
k3_C = r1 * (C + 0.5 * h * k2_C) * (1 - (C + 0.5 * h * k2_C) / K1 - gama1 * W / K2) - ...
       alpha1 * (C + 0.5 * h * k2_C) * I + beta1 * (C + 0.5 * h * k2_C) + sigema4 * B * (C + 0.5 * h * k2_C) + ...
       sigema5 * F * (C + 0.5 * h * k2_C);
k4_C = r1 * (C + h * k3_C) * (1 - (C + h * k3_C) / K1 - gama1 * W / K2) - ...
       alpha1 * (C + h * k3_C) * I + beta1 * (C + h * k3_C) + sigema4 * B * (C + h * k3_C) + ...
       sigema5 * F * (C + h * k3_C);
C = C + (h / 6) * (k1_C + 2 * k2_C + 2 * k3_C + k4_C);


    % 2. dW/dt
k1_W = r2 * W * (1 - W / K2 - gama2 * C / K1) - alpha2 * W * I;
k2_W = r2 * (W + 0.5 * h * k1_W) * (1 - (W + 0.5 * h * k1_W) / K2 - gama2 * C / K1) - ...
       alpha2 * (W + 0.5 * h * k1_W) * I;
k3_W = r2 * (W + 0.5 * h * k2_W) * (1 - (W + 0.5 * h * k2_W) / K2 - gama2 * C / K1) - ...
       alpha2 * (W + 0.5 * h * k2_W) * I;
k4_W = r2 * (W + h * k3_W) * (1 - (W + h * k3_W) / K2 - gama2 * C / K1) - ...
       alpha2 * (W + h * k3_W) * I;
W = W + (h / 6) * (k1_W + 2 * k2_W + 2 * k3_W + k4_W);


    % 3. dI/dt
k1_I = r3 * I * (1 - I / K3) - alpha3 * I * (R + B) - beta3 * I + sigema1 * W *I;
k2_I = r3 * (I + 0.5 * h * k1_I) * (1 - (I + 0.5 * h * k1_I) / K3) - alpha3 * (I + 0.5 * h * k1_I) * (R + B) - ...
       beta3 * (I + 0.5 * h * k1_I) + sigema1 * W * (I + 0.5 * h * k1_I);
k3_I = r3 * (I + 0.5 * h * k2_I) * (1 - (I + 0.5 * h * k2_I) / K3) - alpha3 * (I + 0.5 * h * k2_I) * (R + B) - ...
       beta3 * (I + 0.5 * h * k2_I) + sigema1 * W * (I + 0.5 * h * k2_I);
k4_I = r3 * (I + h * k3_I) * (1 - (I + h * k3_I) / K3) - alpha3 * (I + h * k3_I) * (R + B) - ...
       beta3 * (I + h * k3_I) + sigema1 * W * (I + h * k3_I);
I = I + (h / 6) * (k1_I + 2 * k2_I + 2 * k3_I + k4_I);


    % 4. dB/dt
k1_B = r4 * B * (1 - B / K4 - gama3 * R / K5) + sigema2 * I * B;
k2_B = r4 * (B + 0.5 * h * k1_B) * (1 - (B + 0.5 * h * k1_B) / K4 - gama3 * R / K5) + ...
       sigema2 * I * (B + 0.5 * h * k1_B);
k3_B = r4 * (B + 0.5 * h * k2_B) * (1 - (B + 0.5 * h * k2_B) / K4 - gama3 * R / K5) + ...
       sigema2 * I * (B + 0.5 * h * k2_B);
k4_B = r4 * (B + h * k3_B) * (1 - (B + h * k3_B) / K4 - gama3 * R / K5) + ...
       sigema2 * I * (B + h * k3_B);
B = B + (h / 6) * (k1_B + 2 * k2_B + 2 * k3_B + k4_B);


    % 5. dR/dt
k1_R = r5 * R * (1 - R / K5 - gama4 * B / K4) + sigema3 * I * R;
k2_R = r5 * (R + 0.5 * h * k1_R) * (1 - (R + 0.5 * h * k1_R) / K5 - gama4 * B / K4) + ...
       sigema3 * I * (R + 0.5 * h * k1_R);
k3_R = r5 * (R + 0.5 * h * k2_R) * (1 - (R + 0.5 * h * k2_R) / K5 - gama4 * B / K4) + ...
       sigema3 * I * (R + 0.5 * h * k2_R);
k4_R = r5 * (R + h * k3_R) * (1 - (R + h * k3_R) / K5 - gama4 * B / K4) + ...
       sigema3 * I * (R + h * k3_R);
R = R + (h / 6) * (k1_R + 2 * k2_R + 2 * k3_R + k4_R);


    % 6. dF/dt
k1_F = r6 * F * (1 - F / K6);
k2_F = r6 * (F + 0.5 * h * k1_F) * (1 - (F + 0.5 * h * k1_F) / K6);
k3_F = r6 * (F + 0.5 * h * k2_F) * (1 - (F + 0.5 * h * k2_F) / K6);
k4_F = r6 * (F + h * k3_F) * (1 - (F + h * k3_F) / K6);
F = F + (h / 6) * (k1_F + 2 * k2_F + 2 * k3_F + k4_F);


    % Updating discrete states of species using Markov processes
    for i = 1:6
        current_state = state(i);  
        transition_matrix = transition_matrices{i}; 
        transition_probs = transition_matrix(current_state, :);  
        new_state = find(rand <= cumsum(transition_probs), 1); 
        state(i) = new_state;  
    end
end

% Mapping changes in species populations
figure;
subplot(2, 3, 1); 
area(C_history, 'FaceColor', '#778f83', 'LineWidth', 2, 'EdgeColor', '#778f83', 'FaceAlpha', 0.5);
title('C','FontSize',12);
xticks(0:100:400); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 2); 
area(W_history, 'FaceColor', '#3b66ae', 'LineWidth', 2, 'EdgeColor', '#3b66ae', 'FaceAlpha', 0.5);
title('W','FontSize',12);
xticks(0:100:400); yticks(0:2.5:15); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 3); 
area(I_history, 'FaceColor', '#cac881', 'LineWidth', 2, 'EdgeColor', '#cac881', 'FaceAlpha', 0.5);
title('I','FontSize',12);
xticks(0:100:400); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 4); 
area(B_history, 'FaceColor', '#b0a8b9', 'LineWidth', 2, 'EdgeColor', '#b0a8b9', 'FaceAlpha', 0.5);
title('B','FontSize',12);
xticks(0:100:400); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 5);
area(R_history, 'FaceColor', '#ffc75f', 'LineWidth', 2, 'EdgeColor', '#ffc75f', 'FaceAlpha', 0.5);
title('R','FontSize',12);
xticks(0:100:400); yticks(0:0.5:3); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 6); 
area(F_history, 'FaceColor', '#d65db1', 'LineWidth', 2, 'EdgeColor', '#d65db1', 'FaceAlpha', 0.5);
title('F','FontSize',12);
xticks(0:100:400); 
xlabel('Time   /Day','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Number  / 10 m^2','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;


% Select data points to be plotted at 20 time-step intervals
step_interval = 20;  

% Status (1:Low, 2:Medium, 3:High)
figure;
subplot(2, 3, 1); 
plot(states_history(1:step_interval:end, 1), 'LineWidth', 1, 'Color', '#778f83', 'Marker', 'o', 'MarkerFaceColor', '#778f83');
title('C','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 2); 
plot(states_history(1:step_interval:end, 2), 'LineWidth', 1, 'Color', '#3b66ae', 'Marker', 'o', 'MarkerFaceColor', '#3b66ae');
title('W','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 3); 
plot(states_history(1:step_interval:end, 3), 'LineWidth', 1, 'Color', '#cac881', 'Marker', 'o', 'MarkerFaceColor', '#cac881');
title('I','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 4);
plot(states_history(1:step_interval:end, 4), 'LineWidth', 1, 'Color', '#b0a8b9', 'Marker', 'o', 'MarkerFaceColor', '#b0a8b9');
title('B','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 5); 
plot(states_history(1:step_interval:end, 5), 'LineWidth', 1, 'Color', '#ffc75f', 'Marker', 'o', 'MarkerFaceColor', '#ffc75f');
title('R','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;

subplot(2, 3, 6); 
plot(states_history(1:step_interval:end, 6), 'LineWidth', 1, 'Color', '#d65db1', 'Marker', 'o', 'MarkerFaceColor', '#d65db1');
title('F','FontSize',12);
xlabel('Time   step','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
ylabel('Status','FontName','Times New Roman','Linewidth', 2, 'FontWeight', 'bold','FontSize',12);
grid on;
