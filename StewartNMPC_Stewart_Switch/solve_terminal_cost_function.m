function    P = solve_terminal_cost_function(lambda,MPC)
%MPC model
%lambda > 0 denotes  终端代价函数的系数，如果要求终端状态约束集大一些，而对预测控制消耗的全局性能指标并不是很关心
%可以取较大的lambda，反之，取较小的lambda.默认值取1即可
A = MPC.A;
B = MPC.B;
Q = MPC.Q;
R = MPC.R;
I = eye(size(A,1),size(A,2));
%第二步假设第一步在原点近似求得的线性化系统是稳定的，使用线性二次型调节器（lqr）求得状态线性反馈的最优控制规律
K = lqr(A,B,Q,R);
K = -K;
% K = -[2.0 2.0]
%第三步，构造Lypunuvo函数，求解终端代价函数
Ak = A + B*K;
Qk = Q + K'*R*K;

P = lyap((Ak+lambda*I)',Qk);
%P = P;








