%%
% stewart platform 的几何参数
% length in m
% angles in radians
% mass in kg
% time in seconds

% 底座和平台的半径
rad_upper = 1.0;
rad_lower = 2.0;

% 底座和平台与支腿连接点的角度
ang          = 30;
angles_upper = [ 60-ang, 60-ang, 180-ang, 180-ang, 300-ang, 300-ang] * pi/180;
angles_lower = [ 30-ang, 90-ang, 150-ang, 210-ang, 270-ang, 330-ang ] * pi/180.0;

% zero position
L0 = 464.5435664601383;

% 定义建立（底座与支腿）连接点坐标系的中间单位向量（即：课本上的连接点的坐标单位化ai/|ai|,也是本程序中的rjj/|rjj|）
%%{
U11 = zeros(3, 6);
u10 = [ 0 -1 0 ]';
U11(:,1) = rotationz( angles_lower( 1 ) ) * u10;
U11(:,2) = rotationz( angles_lower( 2 ) ) * u10;
U11(:,3) = rotationz( angles_lower( 3 ) ) * u10;
U11(:,4) = rotationz( angles_lower( 4 ) ) * u10;
U11(:,5) = rotationz( angles_lower( 5 ) ) * u10;
U11(:,6) = rotationz( angles_lower( 6 ) ) * u10;
%%}
% 定义底座与支腿连接点坐标矩阵
raa = zeros(3,6);
r0  = [ 0 -rad_lower 0 ]';
raa(:,1) = rotationz( angles_lower( 1 ) ) * r0;
raa(:,2) = rotationz( angles_lower( 2 ) ) * r0;
raa(:,3) = rotationz( angles_lower( 3 ) ) * r0;
raa(:,4) = rotationz( angles_lower( 4 ) ) * r0;
raa(:,5) = rotationz( angles_lower( 5 ) ) * r0;
raa(:,6) = rotationz( angles_lower( 6 ) ) * r0;

% 定义平台与支腿连接点坐标矩阵
rbb = zeros(3,6);
r0  = [ 0 -rad_upper 0 ]';
rbb(:,1) = rotationz( angles_upper( 1 ) ) * r0;
rbb(:,2) = rotationz( angles_upper( 2 ) ) * r0;
rbb(:,3) = rotationz( angles_upper( 3 ) ) * r0;
rbb(:,4) = rotationz( angles_upper( 4 ) ) * r0;
rbb(:,5) = rotationz( angles_upper( 5 ) ) * r0;
rbb(:,6) = rotationz( angles_upper( 6 ) ) * r0;

% 平台的质量和转动惯量矩阵
mp = 1150;
Ip  = diag( [570, 285, 285]' );

% 下肢（气缸）质量、长度和转动惯量
m1i = 85; 
I1i = diag( [16, 16, 0 ]' );
L1_center = 1.5/2;

% 上肢（活塞）质量、长度和转动惯量
m2i = 22;
I2i = diag( [ 4.1, 4.1, 0 ]' );
L2_center = 1.5/2;

% gravity
gravity = [ 0; 0; -9.8 ];