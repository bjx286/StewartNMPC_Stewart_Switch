function xnew = model_stewart(x, u)


%% Stewart Nonlinear model update
% Dynamic matrix update

%syms xp1 xp2 xp3 rp1 rp2 rp3 vp1 vp2 vp3 wp1 wp2 wp3 Jp_p  M_p  C_p  G_p

% [ Jp_p, M_p, C_p, G_p ] = stewart_dynamic_matrix(x(1:6),x(7:12))

% [ Jp_p, M_p, C_p, G_p ] = stewart_dynamic_matrix(x(1:6),x(7:12));
% [ Jp_p, M_p, C_p, G_p ] = stewart_dynamic_matrix([0 0 3 0 0 0]',[0 0 0 0 0 0]')
% Dynamic equation differential form
    %  X_p = x(1:6);
    %  dX_p = x(7:12);
    %  U_p = M_p*ddX_p + C_p*dX_p + G_p;
    [A_p,B_p,C_p,D_p] = linmod('StewartPlatformPlant',x,[]);
    [A_p,B_p,C_p,D_p] = minreal(A_p,B_p,C_p,D_p); % Obtain a minimal realization of the linearized model.
 
% Dynamic equation state space form
%  A_p = [zeros(6,6) eye(6,6);zeros(6,6)  -M_p'*C_p];
%  B_p = [zeros(6,6);inv(M_p)];
 X = x;
 U_p = [u(1);u(2);u(3);u(4);u(5);u(6)];
 dX = A_p*X + B_p*U_p;
 
 dt = 0.1;                   % sampling time
 xnew = x(:) + dt*dX(:);      % state update