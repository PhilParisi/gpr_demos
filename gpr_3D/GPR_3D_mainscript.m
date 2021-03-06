%%% 3D GPR in MATLAB // URI Phillip Parisi - Update June 2022

% THIS VERSION USES 3D DATA (technically, 2.5D data)
% X and Y --> position on seafloor (training inputs)
% Z --> depth of seafloor at given position (training outputs)

%%%% GUIDE TO USE
%%% .m files you need:
% this script, which is the main script
% gpr_functions folder (located one directory above mainscript)

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
% - Kernel Hyperparameters
% - not recommended to change the map_size, but possible

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP
clc, clearvars, close all, format compact

% Kernel Hyperparameters
hp.L = 3.0;               % lengthscale (high = smoother, low = noisier)
hp.sigma_p = 4.0;         % process noise (output scale / vertical scale)
hp.sigma_n = 0.2;         % sensor noise (used to create W)
hp.kerneltype = 'exact';  % 'exact' or 'sparse' approximate kernel

%%%%%% Generate Random Training Data + Noise (2.5D)
map_size = 20; % side of square map dimensions, keep this even if you can 
x = (-map_size/2:(map_size/2-1));
y = (-map_size/2:(map_size/2-1));
[X,Y] = meshgrid(x,y);

% Matrices map_size x map_size
X = X + randn(map_size,map_size);
Y = Y + randn(map_size,map_size);
Z = 1*sin(X) + Y/20;                    % shape of curve!
% gaussian nosie
noise.mu = 0; noise.sigma = hp.sigma_n;
Z = Z + normrnd(noise.mu,noise.sigma,map_size,map_size)/5;

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
scatter3(X(:,1),X(:,2),Y,'k.')
xlabel('X'),ylabel('Y'),zlabel('Depth')
title('Raw Training Data'), zlim([-5 5])

%%%%%% Prediction Points STAR
% predict at more points than the raw data
predict_pts = map_size*2;
x_star = (-predict_pts/2:0.5:(predict_pts/2)-0.5); 
y_star = (-predict_pts/2:0.5:(predict_pts/2)-0.5); %adding pts twice as dense
[x_star,y_star] = meshgrid(x_star,y_star);
x_star = reshape(x_star,(2*predict_pts)^2,1);
y_star = reshape(y_star,(2*predict_pts)^2,1);

X_Star = [x_star y_star]; % m x 2 % vertical array, can add training data for prediction X
clear x_star y_star

disp('...test data generated, run the next section...')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX CALCS
clc, close all, disp('...running GPR, may take a minute...')
tic 

% Calculate V and Inv(V)                            % depends on training x-points only
W = (hp.sigma_n^2)*eye(nnum);                       % Whitenoise (identity * sigmasquared)
V = K_Function(X,X,hp) + W;                         % Calculate Covariance Matrix using Kernel

% Generate K Parameters
tic
K_Star = K_Function(X_Star,X,hp);                   % Calculate K_Star for New Point(s)
K_StarStar = K_Function(X_Star,X_Star,hp);          % Calculate K_StarStar for New Point(s)
prediction_kernel_calc_time = toc;

% Cholesky Decomposition
L = chol(V,'lower');                                % Lower triangular cholesky factor

% Calculate Predictions!                                    % Finally bring in the training y-points here
Y_Star_Hat = K_Star * CholeskySolve(L,Y);                   % Mean Predictions (mean values of Gaussians)
CapSigma_Star = K_StarStar-K_Star*CholeskySolve(L,K_Star'); % Variance Predictions (prediction covariance matrix)
Y_Star_Var = diag(CapSigma_Star);                           % The diagonals store the variances we want!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOG MARGINAL LIKELIHOOD

% How good is our fit? Use this to tune hyperparameters
LML = calcLML(L,Y,nnum);
AlgoTime = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS & OUTPUTS

%%%% Output LML
clc
fprintf('KernelPredictionTime = %1.2f.\n',prediction_kernel_calc_time)
fprintf('AlgoTime = %1.2f.\n',AlgoTime)
fprintf('Log Marginal Likelihood is %1.1f. Tune hyperparams for better fit.\n',LML)

%%%% Plot Raw Data
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
xlabel('X'),ylabel('Y'),zlabel('Depth')
title('Raw Training Data'), zlim([-5 5])

%%%%% Plot Prediction Data
% 2sigma are plotted rather than variance
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    2, ... % marker size
    'k') % color of pts
scatter3(X_Star(:,1),X_Star(:,2),Y_Star_Hat,...  % Predictions
    10,... % marker size
    2*sqrt(Y_Star_Var), ... % 2sigma
    'filled') % color of pts
xlabel('X'), ylabel('Y'), zlabel('Depth')
title('GPR 3D - Seafloor Ripples')
zlim([-5 5])
legend('Predictions','Location','North')
colormap(); 
bar = colorbar();
ylabel(bar,'variance')
xlabel('X'), ylabel('Y'), zlabel('Depth')
title('GPR 3D - Seafloor Ripples')
zlim([-5 5])

% Training and Prediction Data
figure
scatter3(X(:,1),X(:,2),Y, ...       % Training Data
    5, ... % marker size
    'r') % color of pts
hold on

scatter3(X_Star(:,1),X_Star(:,2),Y_Star_Hat,...  % Predictions
    10,... % marker size
    2*sqrt(Y_Star_Var), ... % 2sigma
    'filled') % color of pts
xlabel('X'), ylabel('Y'), zlabel('Depth')
title('GPR 3D - Seafloor Ripples')
zlim([-5 5])
colormap(); 
legend('Raw Data','Predictions','Location','North')
bar = colorbar();
ylabel(bar,'variance')

