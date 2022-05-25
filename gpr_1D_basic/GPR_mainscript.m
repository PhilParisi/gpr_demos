%%% GPR in MATLAB // URI Phillip Parisi - Update April 2022
%%% Approach from Dr. Kristopher Krasnosky

%%% Version Control Notes
% Added Randomness (process noise)
% Training Data X-Values are Randomized
% Everything Matches KK's Disseration Line 22

%%%% GUIDE TO USE
% .m files you need (all in one folder):
% 1 CalcKernel.m
% 2 CalcKStar.m
% 3 K_Function.m
% 4 this script, which is the main script

%%%% RUNNING & PARAMETERS TO TWEAK
% this script should be good to run out-of-the-box
% there will be randomly generated x-data (training data), plotted in blue
% the GPR data points are plotted in red w/ error bars (uncertainty)

% You can TUNE
% - Kernel Hyperparameters: lengthscale -> hp.L, verticalscale, hp.sig
% - the Training Data + Noise section (nnum is the range for the raw trainig data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP
tic, clc, clearvars, close all, format compact

% Kernel Hyperparameters
hp.L = 4;           % lengthscale (high = smoother, low = noisier)
hp.sig = 3;         % output scale / vertical scale  

% Training Data + Noise (aka Raw Data)
nnum = 40; X_beg = -nnum; X_end = nnum;
X = (X_end - X_beg)*rand(nnum,1) + X_beg;           % vertical array, training X, uniform random
Y = 3*sin(2*pi/40*X) + (rand(nnum,1)*1 - 0.5);      % vertical array, training Y, sinusoidal + noise

% Prediction Points STAR
X_Star = [[(-15+X_beg):2:(15+X_end)]'; X];            % vertical array, add training data for prediction X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS

% Calculate V and Inv(V)                            % depends on training x-points only
var_sig = 0.8;                                      % Process Noise (noise from sensor)
W = (var_sig.^2)*eye(nnum);                         % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                            % Calculate Covariance Matrix using Kernel
V_inv = inv(V);


% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);                      % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);             % Calculate K_StarStar for New Point(s)


% Calculate Predictions!                            % Finally bring in the training y-points here
Y_Star_Hat = K_Star * V_inv * Y;                    % Y Predictions (mean values of Gaussians)
Sigma_Star = K_StarStar - K_Star * V_inv * K_Star'; % Variance Predictions (gives us mean, var for each pt)
Sigma_Star = diag(Sigma_Star);                      % The diagonals store the variances we want!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS

% Plot Training and Prediction Data w/ Error Bars
figure
plot(X_Star,Y_Star_Hat,'ro','MarkerFaceColor','r','MarkerEdgeColor','k'), hold on
errorbar(X_Star,Y_Star_Hat,Sigma_Star,'r.','LineWidth',2)
plot(X,Y,'bo','MarkerFaceColor','b','MarkerSize',8), grid on, hold on
xlabel('X Values'), ylabel('Y Values'), title('Gaussian Process Regression')
legend('Predicted', 'Uncertainty','Raw Data')

toc
