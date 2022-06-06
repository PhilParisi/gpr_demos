%%% 2.5D GPR in MATLAB // URI Phillip Parisi - Update June 2022
%%% Approach from Dr. Kristopher Krasnosky

%%% Version Control Notes
% THIS VERSION USES 3D DATA (technically, 2.5D data)
% X and Y position on seafloor (training inputs)
% Z as depth of seafloor (training outputs,predictions)

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP
clc, clearvars, close all, format compact

% Kernel Hyperparameters
hp.L = 2;           % lengthscale (high = smoother, low = noisier)
hp.sigma = 3;         % output scale / vertical scale  

%%%%%% Generate Random Training Data + Noise (2.5D)
map_size = 20; % side of square map dimensions, keep this even if you can 
x = (-map_size/2:(map_size/2-1));
y = (-map_size/2:(map_size/2-1));
[X,Y] = meshgrid(x,y);

% Matrices map_size x map_size
X = X + randn(map_size,map_size);
Y = Y + randn(map_size,map_size);
Z = 1*sin(X) + Y/20;                    % shape of curve!
Z = Z + randn(map_size,map_size)/5;

% Reshape (map_size^2 x 1)
X1 = reshape(X,map_size^2,1);
Y1 = reshape(Y,map_size^2,1);
Z1 = reshape(Z,map_size^2,1);

% Rename and Clear Variables
X = [X1 Y1];
Y = Z1;
nnum = length(Y);
clearvars -except X Y map_size nnum hp

% Plot
figure
scatter3(X(:,1),X(:,2),Y,'.')
xlabel('X'),ylabel('Y'),zlabel('Depth')
title('Raw Training Data'), zlim([-5 5])

%%%%%% Prediction Points STAR
x_star = (-map_size/2:(map_size/2-1)); 
y_star = (-map_size/2:(map_size/2-1)); 
[x_star,y_star] = meshgrid(x_star,y_star);
x_star = reshape(x_star,map_size^2,1);
y_star = reshape(y_star,map_size^2,1);

X_Star = [x_star y_star]; % m x 2 % vertical array, can add training data for prediction X
clear x_star y_star

disp('...test data generated, run the next section...')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS
clc, close all
tic 

% Calculate V and Inv(V)                            % depends on training x-points only
var_sig = 0.8;                                      % Process Noise (noise from sensor)
W = (var_sig.^2)*eye(nnum);                         % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                         % Calculate Covariance Matrix using Kernel
V_inv = inv(V); %expensive!


% Generate K Parameters
K_Star = K_Function(X_Star,X,hp);                      % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);             % Calculate K_StarStar for New Point(s)


% Calculate Predictions!                            % Finally bring in the training y-points here
Y_Star_Hat = K_Star * V_inv * Y;                    % Y Predictions (mean values of Gaussians)
CapSigma_Star = K_StarStar - K_Star * V_inv * K_Star'; % Variance Predictions (gives us mean, var for each pt)
Y_Star_Var = diag(CapSigma_Star);                      % The diagonals store the variances we want!

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS

% Plot Original Data
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
xlabel('X'),ylabel('Y'),zlabel('Depth')
title('Raw Training Data'), zlim([-5 5])

% Plot Prediction Data
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
scatter3(X_Star(:,1),X_Star(:,2),Y_Star_Hat,...  % Predictions
    10,... % marker size
    Y_Star_Var, ...
    'filled') % color of pts
colormap(); 
bar = colorbar();
ylabel(bar,'variance')
xlabel('X'), ylabel('Y'), zlabel('Depth')
title('GPR 2.5D - Seafloor Ripples')
zlim([-5 5])

% Training and Prediction Data
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
hold on

% Plot Prediction Data
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
scatter3(X_Star(:,1),X_Star(:,2),Y_Star_Hat,...  % Predictions
    10,... % marker size
    Y_Star_Var, ...
    'filled') % color of pts
colormap(); 
bar = colorbar();
ylabel(bar,'variance')
xlabel('X'), ylabel('Y'), zlabel('Depth')
title('GPR 2.5D - Seafloor Ripples')
zlim([-5 5])

