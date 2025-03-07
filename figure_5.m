% figure_5
clc,clear;
% parameter setting
r1 = 1.3;   
K1 = 50;  
C0 = 25;  

% Define the time scale
tspan = 1:1:12;  

% Define the differential equation
logisticODE = @(t, C) r1 * C * (1 - C / K1);

% Solving differential equations using ode45
[t, C] = ode45(logisticODE, tspan, C0);

% Data normalization
data = [C];
data_normalized = data./ repmat(sum(data, 2), 1, size(data, 2));

% Plotting results
m = bar(t, data_normalized );
set(m,'EdgeColor','none','facecolor',[0.53, 0.63, 0.58],'FaceAlpha',.5);

xlabel('Time   /Week','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylabel('Percentage   of   Crop','FontName','Times New Roman','Linewidth', 2,'FontSize',12);

% Modification of the vertical scale to percentage form
yticks(0:0.2:1);
yticklabels(strcat(num2str((0:0.2:1)*100), '%'));
set(gca,'YTickLabel',{'0%','20%','40%','60%','80%','100%',''});

% Adjust the y-axis range to leave a blank space at the top
ylim([0, 1.2]); 

legend('$Crop$','Interpreter',"latex",'Linewidth', 2,'location','north','orientation','horizontal','FontSize',12); 
legend('boxoff')
hold off;