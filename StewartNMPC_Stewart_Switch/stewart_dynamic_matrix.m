function [ Jp, M, C, g ] = stewart_dynamic_matrix( Xp, dXp )

% 此函数计算stewart并联机械臂的动态矩阵（运动方程是根据虚功原理求得的）。
% 
% 输入向量参数：
% 1) Xp = [ xp; rp ]; X是6x1维向量，由平台位置向量xp（3x1）和姿态向量rp（3x1）组成
%
% 2) dXp = [ vp; wp ];Xp是6x1维向量，由平台的线速度向量vp和角速度向量wp组成
%
% 输出矩阵参数：
%
% 1) Jp - a 6x6 matrix.
% Jp 是机械臂的jacobian矩阵.
%
% 2) M - a 6x6 matrix.
% M 并联机械臂的质量矩阵.
%
% 3) C - a 6x6 matrix.
% C 是并联机械臂的科氏力和离心力矩阵coriollis matrix.
%
% 4) g - a 6x1 vector.
% g 是并联机械臂的中立向量gravity vector.
%

%定义叉乘矩阵函数
    function [cmm] = skew(vector)
        cmm = [0 -vector(3) vector(2);...
            vector(3) 0 -vector(1);...
            -vector(2) vector(1) 0];
    end
%定义旋转矩阵
    function [ rotZ ] = rotationz(ang)
        rotZ= [ cos(ang) sin(ang) 0;...
                -sin(ang) cos(ang) 0;...
                      0 0 1];     
    end

    function [ rotY ] = rotationy(ang)
        rotY= [cos(ang) 0 -sin(ang);...
                     0 1 0;...
               sin(ang) 0 cos(ang)];       
    end


    function [ rotX ] = rotationx(ang)
        rotX= [      1 0 0;...
               0 cos(ang) sin(ang);...
               0 -sin(ang) cos(ang)];      
    end


% 检查输入参数
if max( size(Xp) ~= [ 6 1 ] )
    error('Xp must be a 6x1 vector');
end
if max( size(dXp) ~= [ 6 1 ] )
    error('dXp must be a 6x1 vector');
end

% 定义向量
thetap = zeros(3,1);
%rbb    = zeros(3,6);%底座与支腿的连接点坐标
%raa    = zeros(3,6);%平台与支腿的连接点坐标
Jp     = zeros(6,6);
%U11    = zeros(3,6);%底座与支腿的连接点坐标系的中间向量矩阵
Ri     = zeros(3,3);

% 加载并联机械臂的参数数据
%stewart_init_data;
%%
% 输入数据 
r  = Xp(1:3,1);%平台的位置向量
q  = Xp(4:6,1);%平台的姿态（四元数）
R  = angle2dcm(q(1),q(2),q(3),'XYZ');%平台的姿态旋转矩阵
vp = dXp(1:3,1);%平台的线速度
wp = dXp(4:6,1);%平台的角速度

% 初始化输出
M = [ mp*eye(3)  zeros(3,3); zeros(3,3) R*Ip*R' ];

C = [ zeros(3,3) zeros(3,3); zeros(3,3) skew(wp)*R*Ip*R'];

g = [ mp*gravity; zeros(3,1) ];

% 对每只腿
for k = 1:6
    %%%% 逆运动学: 位置
    ri  = R * rbb(:,k);%课本中的bi
    di  = r + ri - raa(:,k) ;%支腿的方向向量
    u1i = U11(:,k);%底座与支腿的连接点坐标系的中间向量
    u3i = di / norm( di ) ;%支腿的单位方向向量
    u2i = skew( u3i ) * u1i; u2i = u2i / norm( u2i );%u2i是与u3i和u1i垂直的单位向量
    u4i = skew( u2i ) * u3i; u4i = u4i / norm( u4i );%...
    c1i = L1_center * u3i;
    c2i = di - L2_center * u3i;
    % 肢体的旋转矩阵
    Ri(:,1) = u4i;
    Ri(:,2) = u2i;
    Ri(:,3) = u3i;
    %%%% 逆运动学: 速度
    rpi       =  skew( wp ) * ri;
    vbi       =  vp + rpi;
    fi        =  u1i' * skew( di )* u2i;
    thetap(1) =  u2i' * vbi / fi;
    thetap(2) =  u1i' * ( u3i*u3i' - eye(3) ) * vbi / fi;
    thetap(3) =  u3i' * vbi;
    w1i       = thetap(1)*u1i + thetap(2)*u2i;
    w2i       = w1i;
    % 支腿的 jacobian 矩阵
    Gi  = u1i * u2i';
    Hi  = u2i*u1i'*u3i*u3i';
    K1i = (Gi - Gi' + Hi) / fi;
    K2i = skew(c2i) * K1i - u3i * u3i';
    % first body jacobian
    J1i = [ -skew(c1i) * K1i   skew(c1i) * K1i * skew(ri); ...
             K1i              -K1i * skew(ri) ];
    % second body jacobian
    J2i = [ -K2i   K2i * skew( ri ); ...
             K1i  -K1i * skew( ri ) ];
    % 机械臂的 jacobian 矩阵
    Jp(k,:)   = [ u3i'  -u3i'*skew(ri) ];
    %% 更新 mass and coriollis matrix
    % 支腿的质量矩阵
    M1i = [ m1i*eye(3)  zeros(3,3); zeros(3,3)  Ri*I1i*Ri' ];
    M2i = [ m2i*eye(3)  zeros(3,3); zeros(3,3)  Ri*I2i*Ri' ];
    % 支腿的 coriollis matrix
    C1i = [ zeros(3,3) zeros(3,3); zeros(3,3) skew(w1i)*Ri*I1i*Ri' ];
    C2i = [ zeros(3,3) zeros(3,3); zeros(3,3) skew(w2i)*Ri*I2i*Ri' ];
    % update matrices and vectors
    M = M + J1i' * M1i * J1i + J2i' * M2i * J2i;
    C = C + J1i' * C1i * J1i + J2i' * C2i * J2i;
    g = g + J1i' * [ m1i*gravity; zeros(3,1) ] + J2i' * [ m2i*gravity; zeros(3,1) ];
    %%%% 计算 dot(J)
    up2i = thetap(1) * skew(u1i) * u2i;
    up3i = (thetap(1) * skew(u1i) + thetap(2) * skew(u2i) ) * u3i;
    dpi  = thetap(3) * u3i + norm(di) * up3i;
    cp1i = L1_center * up3i;
    cp2i = dpi - L2_center * up3i;
    fpi  = u1i' * skew(di) * up2i - u1i'*skew(u2i)*dpi;
    Gpi  = u1i * up2i';
    Hpi  = up2i*u1i'*u3i*u3i' + u2i*u1i'*up3i*u3i' + u2i*u1i'*u3i*up3i';
    Kp1i = (Gpi - Gpi' + Hpi)/fi - fpi*( Gi - Gi' + Hi )/(fi*fi);
    % 下肢 jacobian 矩阵
    Jp1i = [ -skew(cp1i)*K1i-skew(c1i)*Kp1i  skew(cp1i)*K1i*skew(ri)+skew(c1i)*Kp1i*skew(ri)+skew(cp1i)*K1i*skew(rpi); ...
              Kp1i   -Kp1i*skew(ri)-K1i*skew(rpi) ];
    % 上肢 jacobian 矩阵
    Kp2i = skew(cp2i)*K1i + skew(c2i)*Kp1i - up3i*u3i' - u3i*up3i';
    Jp2i = [ -Kp2i  Kp2i*skew(ri)+K2i*skew(rpi);...
             Kp1i  -Kp1i*skew(ri)-K1i*skew(rpi) ];
    %% 更新 corriollis matrix
    C = C + J1i'*M1i*Jp1i + J2i'*M2i*Jp2i;
end

% update sign
g = -g;

return
end

