%AHP
clc,clear;
disp('Please enter the criterion/program level judgment matrix A (nth order).');
A=input('A=');
[n,n]=size(A); 
[V,D]=eig(A); 
% Find the largest eigenvalue and its corresponding eigenvector
tempNum=D(1,1); 
pos=1; 
for h=1:n
    if D(h,h)>tempNum
        tempNum=D(h,h);
        pos=h; 
    end
end
w=abs(V(:,pos)); % Find the eigenvector corresponding to the largest eigenvalue
w=w/sum(w); % normalization process
t=D(pos,pos); % t refers to the largest eigenvalue
disp('Criterion/program layer feature vector (i.e., criterion/program layer weights) w=');disp(w);disp('Criteria layer/program layer Maximum characteristic root t=');disp(t);
% consistency test
CI=(t-n)/(n-1);RI=[0 0 0.52 0.89 1.12 1.26 1.36 1.41 1.46 1.49 1.52 1.54 1.56 1.58 1.59 1.60 1.61 1.615 1.62 1.63];
CR=CI/RI(n);
if CR<0.10
    disp('The consistency of this matrix is acceptable!');
    disp('CI=');disp(CI);
    disp('CR=');disp(CR);
else disp('Consistency validation failed for this matrix, please re-grade!');
end

%The guideline layer can be entered：A = [1, 3, 3; 1/3, 1, 1; 1/3, 1, 1];
%Program layers can be entered：B1 = [1, 1, 3, 2, 7; 1, 1, 3, 2, 7; 1/3, 1/3, 1, 1/4, 5; 1/2, 1/2, 4, 1, 5;1/7, 1/7, 1/5, 1/5, 1];
%Program layers can be entered：B2 = [1, 1/5, 1/3, 2, 1/7;5, 1, 2, 5, 1;3, 1/2, 1, 5, 1/2;1/2, 1/5, 1/5, 1, 1/5;7, 1, 2, 5, 1];
%Program layers can be entered：B3 = [1, 3, 1/3, 1, 2;1/3, 1, 1/5, 1/3, 1;3, 5, 1, 4, 5;1, 3, 1/4, 1, 2;1/2, 1, 1/5, 1/2, 1];

%Total ranking weights of scheme i = scheme layer weights of scheme i * corresponding guideline layer weights and summed to obtain, respectively, the data into Excel to solve the problem
