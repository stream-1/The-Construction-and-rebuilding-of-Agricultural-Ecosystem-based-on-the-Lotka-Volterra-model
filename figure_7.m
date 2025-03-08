% figure_7
clc,clear;
% parameter setting
r2 = 1.6; 
K2 = 50; 
alpha2 = 0.06; 

r3 = 1.2;
K3 = 25; 
alpha3 = 0.1; 

r4 = 0.08; 
K4 = 3; 
gama3 = 0.11; 

r5 = 0.1; 
K5 = 2; 
gama4 = 0.09; 

% Define the time scale
tspan = 1:1:12; 

% initial conditions
W0 = 7; 
I0 = 10; 
B0 = 4; 
R0 = 2; 
y0 = [W0; I0; B0; R0]; 

% Define the differential equation
dydt = @(t, y) [
    r2 * y(1) * (1 - y(1)/K2) - alpha2 * y(1) * y(2);
    r3 * y(2) * (1 - y(2)/K3) - alpha3 * y(2) * (y(3) + y(4));
    r4 * y(3) * (1 - y(3)/K4 - gama3 * y(4)/K5);
    r5 * y(4) * (1 - y(4)/K5 - gama4 * y(3)/K4)
];

% Solving differential equations using ode45
[t, y] = ode45(dydt, tspan, y0);

% Extraction results
W = y(:, 1);
I = y(:, 2);
B = y(:, 3);
R = y(:, 4);

% Data normalization
data = [W, I, B, R];
data_normalized = data./ repmat(sum(data, 2), 1, size(data, 2));

% Plotting results
m = bar(t, data_normalized ,'stacked');
set(m(1),'EdgeColor','none','facecolor',[0.26, 0.45, 0.77],'FaceAlpha',.5);
set(m(2),'EdgeColor','none','facecolor',[0.89, 0.88, 0.57],'FaceAlpha',.8);
set(m(3),'EdgeColor','none','facecolor',[0.37, 0.26, 0.54],'FaceAlpha',.5);
set(m(4),'EdgeColor','none','facecolor',[0.71, 0.55, 0.26],'FaceAlpha',.5);

xlabel('Time   /Week','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylabel('Percentage   of   species ','FontName','Times New Roman','Linewidth', 2,'FontSize',12);

% Modification of the vertical scale to percentage form
yticks(0:0.2:1);
yticklabels(strcat(num2str((0:0.2:1)*100), '%'));
set(gca,'YTickLabel',{'0%','20%','40%','60%','80%','100%',''});

% Adjust the y-axis range to leave a blank space at the top
ylim([0, 1.2]); 

legend('$Weed$','$Insect$','$Bat$','$Bird$','Interpreter',"latex",'Linewidth', 2,'location','north','orientation','horizontal','FontSize',12); 
legend('boxoff')
hold off;