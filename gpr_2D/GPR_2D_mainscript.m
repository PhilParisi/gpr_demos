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

% Kernel Hyperparameters [not optimized/trained]
hp.L = 4;           % lengthscale (high = smoother, low = noisier)
hp.sigma = 3;       % output scale (aka vertical scale)

% Training Data + Gaussian Noise (aka Raw Data)
nnum = 100; X_beg = -nnum; X_end = nnum;
X = (X_end - X_beg)*rand(nnum,1) + X_beg;           % vertical array, training X, uniform random
noise.mu = 0; noise.sigma = 0.25;
Y = 3*sin(2*pi/(0.5*nnum)*X) + (normrnd(noise.mu,noise.sigma,nnum,1)*0.5);  % vertical array, training Y, sinusoidal + noise

% Prediction Points STAR
X_Star = [[(-15+X_beg):2:(15+X_end)]'; X];            % vertical array, add training data for prediction X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                            % depends on training x-points only
var_sig = 0.8;                                      % Process Noise (noise from sensor)
W = (var_sig.^2)*eye(nnum);                         % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                         % Calculate Covariance Matrix using Kernel
V_inv = inv(V);


% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);                      % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);             % Calculate K_StarStar for New Point(s)


% Calculate Predictions!                            % Finally bring in the training y-points here
Y_Star_Hat = K_Star * V_inv * Y;                    % Y Predictions (mean values of Gaussians)
CapSigma_Star = K_StarStar - K_Star * V_inv * K_Star'; % Variance Predictions (gives us mean, var for each pt)
Y_Star_Var = diag(CapSigma_Star);                      % The diagonals store the variances we want!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS

% Plot Training and Prediction Data w/ Error Bars
%   Error bars are typically represented as 2sigma (variance is sigma^2
figure
errorbar(X_Star,Y_Star_Hat,2*sqrt(Y_Star_Var),'r.','LineWidth',2), hold on
plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',4) 
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
legend('Predictions \mu,\sigma^2','Raw Data'), grid on

toc
