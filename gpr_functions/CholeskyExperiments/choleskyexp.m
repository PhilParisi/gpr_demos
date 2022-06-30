% Trying to figure out the in's and out's of cholesky decomposition
clc, clearvars, close all

hp.L = 10; hp.sigma = 1.5;

% Generate Covariance Matrix for Practice
nnum = 70;
X = sort(100*rand(1,nnum)');
Y = 2*X;
V = K_Function(X,X,hp) + (1.2^2)*eye(nnum); % NEED TO ADD ERROR TO ENSURE PosSemiDef

% K_Star for testing product for Y_Star (predictive means)


% Inverse Traditionally
tic
V_inv = inv(V);
prod1 = V_inv*Y;
invtime = toc


% Inverse w/ Cholesky Decomp
    % cholesky decomp doesn't work if terms in V are zero (can't round or
    % set to zero!?!?!)
tic   
L = chol(V,'lower');
prod2 = CholeskySolve(L,Y);
choltime = toc


function x = CholeskySolve(L,b)
% solve Ax = b, but Lx = b per notation. L is cholesky factor chol(J,'lower')
% essentially, if you have U*V^-1*W
% you should L = chol(V,'lower'); then replace w/ U*CholeskySolve(L,W)
    x = (L')\(L\b);
end