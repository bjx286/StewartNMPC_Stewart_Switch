function F_Xf = objFncTerminalSet(u,MPC)

%用外包集序列逼近最大终端状态集
x = [-1.00 -1.00 2.00 -1.57 -1.57 -1.57 -2.00 -2.00 -2.00 -2.00 -2.00 -2.00]';
xnew = model_stewart(x, u);                 % update state via nonlinear model

qx = x.'* MPC.Q * x + [u(1) u(2) u(3) u(4) u(5) u(6)]* MPC.R * [u(1);u(2);u(3);u(4);u(5);u(6)];       % stage cost

F_Xf = qx + xnew.' * MPC.P * xnew;



