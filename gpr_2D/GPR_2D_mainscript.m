%%% GPR in MATLAB // URI Phillip Parisi - Update June 2022
tic, clc, clearvars, close all, format compact

%%%% GUIDE TO USE
% .m files you need (all in one folder):
% 1. this script, which is the main script
% gpr_functions folder
    % 2. SqExpKernel.m, the kernel function
    % 3. K_Function.m, which builds the covariance matrix

% Add gpr_functions to the path
% you can do this manually with addpath(.../filepath/gpr_functions) 
    % if below code does not work
dir_path = cd;
idcs = strfind(dir_path,'/');
func_dir = dir_path(1:idcs(end));
func_dir = strcat(func_dir,"gpr_functions");
addpath(func_dir);

%%%% RUNNING & PARAMETERS TO TWEAK
% this script should be good to run out-of-the-box
% there will be randomly generated x-data (training data), plotted in blue
% the GPR data points are plotted in red w/ error bars (uncertainty)

% You can TUNE
% - Kernel Hyperparameters: lengthscale -> hp.L, verticalscale, hp.sigma
% - nnum, number of points you want in the dataset

% Add gpr_functions to the path
% you can do this manually with addpath(.../filepath/gpr_functions) 
    % if below code does not work


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

% Kernel Hyperparameters [not optimized/trained] & Noise
hp.L = 10;                  % lengthscale (high = smoother, low = noisier)
hp.sigma_p = 2;           % process noise (aka vertical scale, output scale)
hp.sigma_n = .6;           % sensor noise (used to create W)
hp.kerneltype = 'exact';    % 'exact' or 'sparse' approximate kernel
                 % Note: for sparse, hp.L = 40 and hp.sigma_p = 2 for hp.sigma_n = 0.6
                 % Note: for exact, hp.L = 10 and hp.sigma_p = 2 for hp.sigma_n = 0.6

% Generate Training Data w/ Gaussian Noise (aka Raw Data)
nnum = 100; X_beg = -nnum; X_end = nnum;
X = sort((X_end - X_beg)*rand(nnum,1) + X_beg);           % vertical array, training X, uniform random
noise.mu = 0; noise.sigma = 0.25;
Y = 3*sin(2*pi/(0.5*nnum)*X) + (normrnd(noise.mu,noise.sigma,nnum,1)*0.5);  % vertical array, training Y, sinusoidal + noise

% Prediction Points (X_Star)
X_Star = [[(-15+X_beg):2:(15+X_end)]'; X];            % vertical array, add training data for prediction X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                                   % depends on training x-points only
W = (hp.sigma_n^2)*eye(nnum);                              % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp,hp.kerneltype) + W;                  % Calculate Covariance Matrix using Kernel
%V_inv = inv(V);


% Generate K Parameters
K_Star = K_Function(X_Star,X,hp,hp.kerneltype);            % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp,hp.kerneltype);   % Calculate K_StarStar for New Point(s)

% Cholesky Decomposition
L = chol(V,'lower');                                       % Lower triangular cholesky factor

% Calculate Predictions!                                   % Finally bring in the training y-points here
%Y_Star_Hat = K_Star * V_inv * Y;                          % Y Predictions (mean values of Gaussians)
Y_Star_Hat = K_Star * CholeskySolve(L,Y);
%CapSigma_Star = K_StarStar - K_Star * V_inv * K_Star';    % Variance Predictions (gives us mean, var for each pt)
CapSigma_Star = K_StarStar-K_Star*CholeskySolve(L,K_Star');
Y_Star_Var = diag(CapSigma_Star);                         % The diagonals store the variances we want!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOG MARGINAL LIKELIHOOD
% How good is our fit?
% note log in compsci means natural log, and matlab's log() is natural log

%LML = -0.5*log(det(V)) - 0.5*Y'*V_inv*Y - 0.5*nnum*log(2*pi)
LML = -0.5*log(det(V)) - 0.5*Y'*CholeskySolve(L,Y) - 0.5*nnum*log(2*pi);
LML

%LML_array=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS

% Organize Data to Plot It (created new obj
sortobj = [X_Star, Y_Star_Hat, Y_Star_Var];
sortobj = sortrows(sortobj);
sorted.X_Star = sortobj(:,1); sorted.Y_Star_Hat = sortobj(:,2); sorted.Y_Star_Var = sortobj(:,3);

% Error Bar Plot (Training Data + Predictoins w/ Error Bars)
    %Error bars are typically represented as 2sigma (variance is sigma^2
% figure
% errorbar(sorted.X_Star,sorted.Y_Star_Hat,2*sqrt(sorted.Y_Star_Var),'r.','LineWidth',2), hold on
% plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',4) 
% xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
% legend('Predictions \mu,\sigma^2','Raw Data'), grid on

% Bounded Plot (Training Data + Predictions + 2sigma Upper and Lower Bound
figure
p1 = plot(sorted.X_Star,sorted.Y_Star_Hat + 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); hold on %upper bound
plot(sorted.X_Star,sorted.Y_Star_Hat - 2*sqrt(sorted.Y_Star_Var),'r','LineWidth',2); % lower bound
p2 = plot(sorted.X_Star,sorted.Y_Star_Hat,'r--','Linewidth',2); % prediction means
p3 = plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',4); % training data
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
grid on, legend([p3 p2 p1],"Training Data","Prediction \mu","Prediction 2\sigma")

toc
