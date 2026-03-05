clc; clear; close all;

addpath('godlike');
addpath('Linearization');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of approximate nonlinear MPC for stewart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 初始化

nx = 12;             % 状态变量的数目
nu = 6;             % 控制输出的数目

Ns = 100;           % 状态域X上样点数目

MPC.Np = 15;        % MPC 预测域长度（预测步数）
MPC.Nu = nu;
MPC.Nx = nx;

%第一步：对stewart_model在[0 0 2 0 0 0]处进行线性化处理
% Linearization.
StewartPlatformSetup
sim('StewartPlatformEquilibrium'); nomForces = Forces'; % 求解平衡点处的力nomForces
[A,B,C,D] = linmod('StewartPlatformPlant',[],nomForces); % 线性化
[A,B,C,D] = minreal(A,B,C,D); % 最小线性化模型

% MPC weight parameters
% a21  =  [9.0418 -4.7826 1.4778 -2.8144 -2.0513 2.4682;...
%          -4.7826 9.0418 2.4682 -2.0513 -2.8144 1.4778;...
%          -2.0513 2.4682 9.0418 -4.7826 1.4778 -2.8144;...
%          -2.8144 1.4778 -4.7826 9.0418 2.4682 -2.0513;...
%          1.4778 -2.8144 -2.0513 2.4682 9.0418 -4.7826;...
%          2.4682 -2.0513 -2.8144 1.4778 -4.7826 9.0418...
%          ];
% b21  =  [0.0002 0.0002 -0.0010 0.0002 -0.0001 0.0021;...
%          0.0002 -0.0010 0.0002 0.0002 0.0021 -0.0001;...
%          -0.0010 0.0002 -0.0001 0.0021 0.0002 0.0002;...
%          0.0002 0.0002 0.0021 -0.0001 0.0002 -0.0010;...
%          -0.0001 0.0021 0.0002 0.0002 -0.0010 0.0002;...
%          0.0021 -0.0001 0.0002 -0.0010 0.0002 0.0002];
% c12  =  [0 0 0 0 0 1;...
%          0 0 0 0 1 0;...
%          0 0 0 1 0 0;...
%          0 0 1 0 0 0;...
%          0 1 0 0 0 0;...
%          1 0 0 0 0 0];
% MPC.A = [zeros(6,6) eye(6,6);a21 zeros(6,6)];
% MPC.B = [zeros(6,6);b21];
% MPC.C = [c12 zeros(6,6)];
%MPC.P = [16.5926, 11.5926; 11.5926, 16.5926];
%stewart 模型参数 parameters
MPC.A = A;
MPC.B = B;
MPC.C = C;
MPC.D = D;
MPC.Q = diag([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]);
MPC.R = diag([1.0,1.0,1.0,1.0,1.0,1.0]);
MPC.P = solve_terminal_cost_function(0.95,MPC);


% 控制量边界
MPC.Ulb = 0*zeros(1,nu);
MPC.Uub = 60000*ones(1,nu);

% 状态边界
% MPC.Xlb = [-0.6 -0.6 +2 -1 -1 -1 -0.060 -0.08 -2 -2 -2 -2]';
% MPC.Xub = [+0.6 +0.6 +3 +1 +1 +1 +0.060 +0.08 +2 +2 +2 +2]';

MPC.Xlb = [-0.6 -0.6 +2 -1 -1 -1 -0.05 -0.04 -0.04 -0.05 -0.04 -0.04]';
MPC.Xub = [+0.6 +0.6 +3 +1 +1 +1 +0.05 +0.04 +0.04 +0.05 +0.04 +0.04]';

%% 选择状态样本
H = haltonset(nx);
H = scramble(H, 'RR2');
X = MPC.Xlb + (MPC.Xub - MPC.Xlb).* net(H, Ns)';     % allocating memory for state samples
clear H

%% Computing optimal control samples

% U = zeros(Ns, nu);  % allocating memory for control samples
% F = zeros(Ns, 1);   % allocating memory for feasibility information

% for k = 1:Ns
%     fprintf('\nIteration: %d of %d', k, Ns);
%     [U(k), F(k)] = find_optimal_NMPC(X(k,:), MPC);
% end


%% Saving samples to save time
save('dataset.mat');