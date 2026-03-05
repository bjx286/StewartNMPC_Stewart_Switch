% function [uk, feas] = find_optimal_NMPC(xk, MPC)
function [uk,feas] = find_optimal_NMPC(xk, MPC)

% for k = 1:size(xk,1)
     xk = xk';
    %xk = [0; 0;3;0;0;0;0;0;0;0;0;0];
    % % Solve the finite-horizon optimal control problem
    % u_opt = GODLIKE(@(u) objFncTerminalSet(u, x0, MPC),...
    %                 MPC.Ulb*ones(MPC.Np,1), MPC.Uub*ones(MPC.Np,1));
    % uk = u_opt(1:MPC.Nu);
    % syms u
    % F_Xf = objFncTerminalSet(u, x0, MPC)
%      u = sym('u',[6 1]);
%      
% %     xnew = [x1;x2];
% %     x = [x1;x2;u];
%     xnew = model_stewart(xk, u);                 % update state via nonlinear model
% %     f_x = w_svm'*xnew + b_svm
% %     solve('xnew','u')
%     qx = xk.'* MPC.Q * xk + u.'* MPC.R * u;       % stage cost
%     F_Xf = qx + xnew.' * MPC.P * xnew; 
%     diff_u1 = diff(F_Xf,u(1));
%     diff_u2 = diff(F_Xf,u(2));
%     diff_u3 = diff(F_Xf,u(3));
%     diff_u4 = diff(F_Xf,u(4));
%     diff_u5 = diff(F_Xf,u(5));
%     diff_u6 = diff(F_Xf,u(6));
%     
%     [u(1),u(2),u(3),u(4),u(5),u(6)] = solve(diff_u1 == 0,diff_u2 == 0,diff_u3 == 0,diff_u4 == 0,diff_u5 == 0,diff_u6 == 0,u(1),u(2),u(3),u(4),u(5),u(6))
%     vpa(u,9)
%     fun = @(u)28562.9*u(3) - 25923.5*u(2) - 6922.41*u(1) - 7098.16*u(4) - 58192.1*u(5) + 41453.2*u(6)...
%               - 1.0*(4.10971e-7*u(1) + 2.65956e-7*u(2) - 2.4395e-7*u(3) + 9.65988e-7*u(4) + 0.000287357*u(5) - 3.61038e-7*u(6)+ 29779.2)...
%               *(1438.35*u(1) - 2924.56*u(2) + 1466.59*u(3) + 7718.23*u(4) - 12372.9*u(5) + 1453.78*u(6) - 3.19426e12)...
%               - 1.0*(0.0000785794*u(2) - 9.064e-8*u(1) + 1.63515e-7*u(3) + 4.53698e-7*u(4) + 2.65956e-7*u(5) + 5.04708e-8*u(6) + 14170.3)...
%               *(702.119*u(1) - 3385.27*u(2) + 2514.34*u(3) + 4460.59*u(4) - 10716.9*u(5) + 2951.61*u(6) - 2.77613e12)...
%               - 1.0*(715.403*u(2) - 3375.24*u(1) + 1446.52*u(3) - 9032.0*u(4) + 5275.23*u(5) + 5080.63*u(6) + 1.67213e12)...
%               *(0.0000784823*u(1) - 9.064e-8*u(2) + 9.35656e-8*u(3) - 2.00283e-7*u(4) + 4.10971e-7*u(5) - 1.96869e-7*u(6) + 17709.4)...
%               - 1.0*(4.53698e-7*u(2) - 2.00283e-7*u(1) + 9.37046e-8*u(3) + 0.000243006*u(4) + 9.65988e-7*u(5) - 1.64071e-7*u(6) - 46508.5)...
%               *(1455.73*u(2) - 2908.5*u(1) + 704.554*u(3) - 10436.3*u(4) + 9139.28*u(5) + 2949.31*u(6) + 2.60175e12)...
%               + (1.96869e-7*u(1) - 5.04708e-8*u(2) + 1.45934e-8*u(3) + 1.64071e-7*u(4) + 3.61038e-7*u(5) - 0.000158434*u(6) - 10527.2)...
%               *(2516.11*u(1) + 1460.67*u(2) - 2916.4*u(3) + 4532.49*u(4) + 2689.33*u(5) - 6844.44*u(6) + 4.22122e11)...
%               + u(1)^2 + u(2)^2 + u(3)^2 + u(4)^2 + u(5)^2 + u(6)^2 ...
%               - 1.0*(1465.38*u(1) + 2514.34*u(2) - 3377.37*u(3) + 2253.74*u(4) + 5391.61*u(5) - 5916.01*u(6) + 1.25678e12)...
%               *(9.35656e-8*u(1) + 1.63515e-7*u(2) + 0.0000783956*u(3) + 9.37046e-8*u(4) - 2.4395e-7*u(5) - 1.45934e-8*u(6) - 17120.9) - 1.44156e13
%         x0=[10,0,10,0,10,0];%优化起点
%     A=[];
%     b=[];
%     Aeq=[];
%     beq=[];
%     lb=[-inf,-inf,-inf,-inf,-inf,-inf];
%     ub=[inf,inf,inf,inf,inf,inf];
%     fval = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
%     %%%%test
%     fun = @(u)0.26*u(1)^2 - 0.11*u(1)*u(2) - 0.23*u(1)*u(3) + 1.4*u(1)*u(4) - 0.82*u(1)*u(5) - 0.8*u(1)*u(6) - 2.6e8*u(1) +...
%               0.27*u(2)^2 -0.39*u(2)*u(3) - 0.7*u(2)*u(4) + 1.7*u(2)*u(5) - 0.46*u(2)*u(6) + 4.4e8*u(2) + ...
%               0.26*u(3)^2 - 0.35*u(3)*u(4) - 0.85*u(3)*u(5) + 0.93*u(3)*u(6) - 2.0e8*u(3) + ...
%               2.5*u(4)^2 - 4.4*u(4)*u(5) - 1.4*u(4)*u(6) - 1.3e9*u(4) + ...
%               3.5*u(5)^2 - 0.85*u(5)*u(6) + 1.8e9*u(5) + ...
%               1.1*u(6)^2 - 1.3e8*u(6)...
%               + 2.4e17 ;
%     
%     clc;clear;clear all; 
%         StewartPlatformSetup
%     sim('StewartPlatformEquilibrium'); nomForces = Forces'; % Extract the equilibrating forces. 
%     [A,B,C,D] = linmod('StewartPlatformPlant',[],nomForces); % Linearize the model about the equilibrium configuration.
%     [A,B,C,D] = minreal(A,B,C,D); % Obtain a minimal realization of the linearized model.      
%           
%     MPC.A = A;
%     MPC.B = B;
%     MPC.C = C;
%     MPC.D = D;
%     MPC.Q = diag([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]);
%     MPC.R = diag([1.0,1.0,1.0,1.0,1.0,1.0]);
%     MPC.P = solve_terminal_cost_function(0.95,MPC);
%     
%     
%     xk = [-1.00 -1.00 2.00 -1.57 -1.57 -1.57 -2.00 -2.00 -2.00 -2.00 -2.00 -2.00]';
%     
%     fun = @objFncTerminalSet;
    
    
%     u = sym('u',[6 1]);
%     dd = model_stewart(xk, u)

%%
    [ Jp_p, M_p, C_p, G_p ] = stewart_dynamic_matrix(xk(1:6),xk(7:12));
    % [ Jp_p, M_p, C_p, G_p ] = stewart_dynamic_matrix([0 0 3 0 0 0]',[0 0 0 0 0 0]')
    % Dynamic equation differential form
        %  X_p = x(1:6);
        %  dX_p = x(7:12);
        %  U_p = M_p*ddX_p + C_p*dX_p + G_p;
    % Dynamic equation state space form
     A_p = [zeros(6,6) eye(6,6);zeros(6,6)  -M_p'*C_p];
     B_p = [zeros(6,6);inv(M_p)];
     X = xk;
    %  U_p = [u(1);u(2);u(3);u(4);u(5);u(6)] - G_p;
    %  dX = A_p*X + B_p*U_p;
    %  
    %  dt = 0.1;                   % sampling time
    %  xnew = xk(:) + dt*dX(:);      % state update
%%
    fun = @(u) xk.'* MPC.Q * xk + [u(1) u(2) u(3) u(4) u(5) u(6)]* MPC.R * [u(1);u(2);u(3);u(4);u(5);u(6)] +...
               (xk(:) + 0.1*(A_p*X + B_p*([u(1);u(2);u(3);u(4);u(5);u(6)])))' * MPC.P * (xk(:) + 0.1*(A_p*X + B_p*([u(1);u(2);u(3);u(4);u(5);u(6)]))) ;
    x0=[0,0,0,0,0,0];%优化起点
    if MPC.w_old == zeros(1,MPC.Nx) & MPC.b_old == 0
        Aneq=[];
        bneq=[];
    else
        U = sym('U',[6 1]);
        xnew = model_stewart(xk,U);
        equation = MPC.w_old*xnew + MPC.b_old;
        for k=1:MPC.Nu
            Aneq(1,k) = double(diff(equation,U(k)));%约束为大于等于0，转化为小于等于时，Aneq前加负号
        end
        bneq=double(subs(equation,U,zeros(MPC.Nu,1)));
    end

    Aeq=[];
    beq=[];
    [u_min,fval] = fmincon(fun,x0,Aneq,bneq,Aeq,beq,MPC.Ulb,MPC.Uub);
    uk = u_min;

    Fx = xk.' * MPC.P * xk;
    minF_Xf = fval;
    if Fx >=  minF_Xf
        feas = +1;                % feasible标签
    else
        feas = -1;                  % infeasible标签
    end

% end
            
% Check if the final state is in the terminal region
% for k = 1:MPC.Np
%     xk = model_test(xk, u_opt(k));         % update state via nonlinear model
% end
% 
% if xk.' * MPC.P * xk <= 0.7 
%     feas = +1;                  % label as feasible
% else
%     feas = -1;                  % label as infeasible
% end