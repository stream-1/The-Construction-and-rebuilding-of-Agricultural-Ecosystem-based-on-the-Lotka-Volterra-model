% figure_3
clc,clear;
% parameter setting
r1 = 1.3; 
K1 = 50; 
theta1 = 0.67;
gama1 = 1.3; 
alpha1 = 0.06; 
beta1 = 0.9; 

r2 = 1.6; 
K2 = 25; 
gama2 = 0.2; 
alpha2 = 0.06; 
beta2 = 0.8; 

r3 = 1.2;
K3 = 25; 
beta3 = 1;

% Define the time scale
tspan = 1:1:10; 

% initial conditions
C0 = 25; 
W0 = 11; 
I0 = 16; 
y0 = [C0; W0; I0]; 

% Define the differential equation
dydt = @(t, y) [
    r1 * y(1) * (1 - y(1)/(theta1*K1) - gama1 * y(2)/K2) - alpha1 * y(1) * y(3) + beta1 * y(1);
    r2 * y(2) * (1 - y(2)/K2 - gama2 * y(1)/(theta1*K1)) - alpha2 * y(2) * y(3) - beta2 * y(2);
    r3 * y(3) * (1 - y(3)/K3) - beta3 * y(3)
];

% Solving differential equations using ode45
[t, y] = ode45(dydt, tspan, y0);

% Extraction results
C = y(:, 1);
W = y(:, 2);
I = y(:, 3);

% Plotting results
data = [C, W, I]; 
m = bar(t, data, "grouped");

set(m(1),'EdgeColor','none','facecolor',[0.53, 0.63, 0.58],'FaceAlpha',.5);
set(m(2),'EdgeColor','none','facecolor',[0.26, 0.45, 0.77],'FaceAlpha',.5);
set(m(3),'EdgeColor','none','facecolor',[0.89, 0.88, 0.57],'FaceAlpha',.8);

legend('$Crop$','$Weed$','$Insect$','Interpreter',"latex",'Linewidth', 2,'FontSize',12); 
legend('boxoff')
xlabel('Time   /Year','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylabel('Number   of   species   / 10 m^2 ','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
legend;
hold off;