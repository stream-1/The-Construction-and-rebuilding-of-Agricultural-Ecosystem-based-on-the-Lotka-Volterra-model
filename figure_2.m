%figure_2
% With a spacing of 40 cm between plants and 50 cm between rows, it was calculated that 50 wheat plants were planted in 10 square meters, which is the same as the spacing between plants and rows in North China.
clc,clear;
%parameter setting
r1 = 1.3;   
K1 = 50; 
C0 = 10;  

% Define the time scale
tspan = 1:1:10;  

% Define the differential equation
logisticODE = @(t, C) r1 * C * (1 - C / K1);

% Solving differential equations using ode45
[t, C] = ode45(logisticODE, tspan, C0);

% Plotting results
figure;
x=bar(t, C);
set(x,'EdgeColor','none','facecolor',[0.53, 0.63, 0.58],'FaceAlpha',.5)
xlabel('Time   /Year','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylabel('Number   of   Crop   / 10 m^2 ','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
legend('$Crop$','Interpreter',"latex",'Linewidth', 2,'FontSize',12);
legend('boxoff')