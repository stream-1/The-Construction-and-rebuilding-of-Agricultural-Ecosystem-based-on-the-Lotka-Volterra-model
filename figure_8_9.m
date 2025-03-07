% figure 8 and figure 9
clear,clc;
% parameter setting
r1 = 1.3; 
K1 = 50; 
theta1 = 0.85;
gama1 = 1.3; 
alpha1 = 0.06; 
beta1 = 0.8; 

r2 = 1.6; 
K2 = 25; 
gama2 = 0.2; 
alpha2 = 0.06; 
beta2 = 0.6; 

r3 = 1.2;
K3 = 25; 
alpha3 = 0.1; 
beta3 = 0.5;
sigema1=0.08;

r4 = 0.08; 
K4 = 3; 
gama3 = 0.11; 
sigema2=0.08;

r5 = 0.1; 
K5 = 2; 
gama4 = 0.09; 
sigema3=0.08;

% Define the time scale
tspan = 0:0.2:24; 

% initial conditions
C0 = 25; 
W0 = 11; 
I0 = 16; 
B0 = 2; 
A0 = 1; 
y0 = [C0; W0; I0; B0; A0]; 

%% Define the system of differential equations 1
dydt1 = @(t, y) [
    r1 * y(1) * (1 - y(1)/(theta1*K1) - gama1 * y(2)/K2) - alpha1 * y(1) * y(3) + beta1 * y(1);
    r2 * y(2) * (1 - y(2)/K2 - gama2 * y(1)/(theta1*K1)) - alpha2 * y(2) * y(3) - beta2 * y(2);
    r3 * y(3) * (1 - y(3)/K3) - alpha3 * y(3) * (y(4) + y(5)) - beta3 * y(3)+ sigema1 * y(2) * y(3);
    r4 * y(4) * (1 - y(4)/K4 - gama3 * y(5)/K5)+ sigema2 * y(3) * y(4) ;
    r5 * y(5) * (1 - y(5)/K5 - gama4 * y(4)/K4)+ sigema3 * y(3) * y(5) 
];

% Solving Systems of Differential Equations 1
[t, y1] = ode45(dydt1, tspan, y0);

% Extraction result 1
C1 = y1(:, 1);
W1 = y1(:, 2);
I1 = y1(:, 3);
B1 = y1(:, 4);
R1 = y1(:, 5);

data1 = [C1, W1, I1, B1, R1];

%% Define the system of differential equations 2
dydt2 = @(t, y) [
    r1 * y(1) * (1 - y(1)/(theta1*K1) - gama1 * y(2)/K2) - alpha1 * y(1) * y(3) + beta1 * y(1);
    r2 * y(2) * (1 - y(2)/K2 - gama2 * y(1)/(theta1*K1)) - alpha2 * y(2) * y(3) ;
    r3 * y(3) * (1 - y(3)/K3) - alpha3 * y(3) * (y(4) + y(5)) - beta3 * y(3)+ sigema1 * y(2) * y(3);
    r4 * y(4) * (1 - y(4)/K4 - gama3 * y(5)/K5)+ sigema2 * y(3) * y(4) ;
    r5 * y(5) * (1 - y(5)/K5 - gama4 * y(4)/K4)+ sigema3 * y(3) * y(5) 
];

% Solving Systems of Differential Equations 2
[t, y2] = ode45(dydt2, tspan, y0);

% Extraction result 2
C2 = y2(:, 1);
W2 = y2(:, 2);
I2 = y2(:, 3);
B2 = y2(:, 4);
R2 = y2(:, 5);
data2 = [C2, W2, I2, B2, R2];

% the result of subtraction
data=data2-data1;

%% Plotting results

plot(t,data(:,1),'Linewidth',2,'Color','#778f83');
hold on;
plot(t,data(:,2),'Linewidth',2,'Color','#3b66ae');
plot(t,data(:,3),'Linewidth',2,'Color','#cac881');
plot(t,data(:,4),'Linewidth',2,'Color','#845ec2');
plot(t,data(:,5),'Linewidth',2,'Color','#ffc75f');
 
xlabel('Time   /Week','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylabel('Number   change   in   species   / 10 m^2','FontName','Times New Roman','Linewidth', 2,'FontSize',12);
ylim([-30,20]);
legend('$Crop$','$Weed$','$Insect$','$Bat$','$Bird$','Interpreter',"latex",'Linewidth', 2,'location','north','orientation','horizontal','FontSize',12); 
legend('boxoff')
hold off;
grid on;


%% anoval test and p-value test
% Setting up group labels
group1 = repmat(1, 121, 1);  
group2 = repmat(2, 121, 1);  

% Merge data and group labels
data_combined = [data1; data2];  
group_combined = [group1; group2];  

% Used to store the p-value for each column
p_values = zeros(1, 5);

% Perform anova1 test for each column
for i = 1:5
    % For each column of data, the anova1 test is performed separately.
    [p_values(i), ~, stats] = anova1(data_combined(:, i), group_combined);
    
    % Display the p-value for each column
    disp(['P-value for column ', num2str(i), ': ', num2str(p_values(i))]);
end

% Significant differences exist if the p-value is less than 0.05
