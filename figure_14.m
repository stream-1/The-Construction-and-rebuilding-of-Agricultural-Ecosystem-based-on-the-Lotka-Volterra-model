% figure_14
clc,clear;
% Create a widescreen graphics window, resize it
figure('Position', [100, 100, 1200, 800]); 

% parameter setting
r1 = 1.3; 
K1 = 50; 
theta1=0.85;
gama1 = 1.3; 
alpha1 = 0.06; 
beta1 = 0.8; 

r2 = 1.6; 
K2 = 25;  
gama2 = 0.2; 
alpha2 = 0.06; 
beta2 = 0.6; 

r3_values=[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]; 
K3 = 25; 
alpha3 = 0.1; 
beta3 = 0.5;

r4 = 0.08; 
K4 = 3; 
gama3 = 0.11; 

r5 = 0.1; 
K5 = 2; 
gama4 = 0.09;  

% Define the time scale
tspan = 2:2:24; 

% initial conditions
C0 = 25; 
W0 = 11; 
I0 = 16; 
B0 = 2; 
R0 = 1; 
y0 = [C0; W0; I0; B0; R0]; 


for i = 1:9
    % Modify the pest growth rate constant r3
    r3 = r3_values(i);
    
    % Define a system of differential equations
    dydt = @(t, y) [
        r1 * y(1) * (1 - y(1)/(K1*theta1) - gama1 * y(2)/K2) - alpha1 * y(1) * y(3) + beta1 * y(1);
        r2 * y(2) * (1 - y(2)/K2 - gama2 * y(1)/(K1*theta1)) - alpha2 * y(2) * y(3) - beta2 * y(2);
        r3 * y(3) * (1 - y(3)/K3) - alpha3 * y(3) * (y(4) + y(5)) - beta3 * y(3);
        r4 * y(4) * (1 - y(4)/K4 - gama3 * y(5)/K5);
        r5 * y(5) * (1 - y(5)/K5 - gama4 * y(4)/K4)
    ];
    
    % Solving Systems of Differential Equations
    [t, y] = ode45(dydt, tspan, y0);
    
    % Extraction results
    C = y(:, 1);
    W = y(:, 2);
    I = y(:, 3);
    B = y(:, 4);
    R = y(:, 5);

    % Data normalization
    data = [C, W, I, B, R];
    data_normalized = data ./ repmat(sum(data, 2), 1, size(data, 2));

    % Create subgraphs with 3 rows and 3 columns layout
    subplot(3, 3, i);
    
    % Plotting stacked bar charts
    m = bar(t, data_normalized, 'stacked');
    set(m(1), 'EdgeColor', 'none', 'FaceColor', [0.53, 0.63, 0.58], 'FaceAlpha', 0.5);
    set(m(2), 'EdgeColor', 'none', 'FaceColor', [0.26, 0.45, 0.77], 'FaceAlpha', 0.5);
    set(m(3), 'EdgeColor', 'none', 'FaceColor', [0.89, 0.88, 0.57], 'FaceAlpha', 0.8);
    set(m(4), 'EdgeColor', 'none', 'FaceColor', [0.37, 0.26, 0.54], 'FaceAlpha', 0.5);
    set(m(5), 'EdgeColor', 'none', 'FaceColor', [0.71, 0.55, 0.26], 'FaceAlpha', 0.5);
    
    % Setting up labels and titles
    xlabel('Time / Week', 'FontName', 'Times New Roman', 'LineWidth', 2,'FontSize',12);
    ylabel('Percentage of species', 'FontName', 'Times New Roman', 'LineWidth', 2,'FontSize',12);
    title(['r_3 = ', num2str(r3)], 'FontName', 'Times New Roman','FontSize',12);
     xticklabels(2:2:24);
    
    % Modification of the vertical scale to percentage form
    yticks(0:0.2:1);
    yticklabels(strcat(num2str((0:0.2:1)*100), '%'));
    set(gca, 'YTickLabel', {'0%', '20%', '40%', '60%', '80%', '100%'});
    
end

% Add legends to the outside of the 9 figures
lgd = legend('$Crop$', '$Weed$', '$Insect$', '$Bat$', '$Bird$', 'Interpreter', 'latex', ...
    'LineWidth', 0.5, 'Location', 'northeastoutside', 'Orientation', 'horizontal','FontSize',12);
set(lgd, 'Position', [0.5 - lgd.Position(3)/2, 0.92, lgd.Position(3), lgd.Position(4)]);


hold off;

% Adjust 'PaperPosition' when saving images
set(gcf, 'PaperPositionMode', 'auto');  
set(gcf, 'PaperSize', [12 8]);  
