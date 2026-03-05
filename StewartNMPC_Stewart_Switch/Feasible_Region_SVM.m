clc; clear; close all;
%load dataset
addpath('godlike');
initFnc;
X = X';
Xfk = X;%此处X为行向量矩阵，以行为单位
MPC.w_old = zeros(1,nx);
MPC.b_old = 0;
h = 0.02; % 网格步长
[X1,X2] = meshgrid(-1:h:1, -1:h:1);
dd = ones(101*101,1);
Ns = size(X,1);
F = ones(Ns, 1);
U = zeros(Ns,nu);
xfi_num = 0;
epslion = 0.01;%循环结束的阈值条件参数，这里取1，其值越小，逼近精度越高，循环次数也越多
Stewart_init_data;
while true
    %Xfk为svm数据训练集，F为训练集标签
     for k = 1:Ns
        fprintf('\nThe %d Iteration: %d of %d',xfi_num, k, Ns);
        %先判断x（k，：）是否在Xfi里面
        if F(k) == 1
            if MPC.w_old * X(k,:)' + MPC.b_old < 0
                F(k) = -1;
            else %如果x（k，：）在Xfi里面，则判断xnew在Xfi里面的约束边界条件，w*xnew + b > 0 和 MPC.Ulb < u < MPC.Uub
                [U(k,:),F(k)] = find_optimal_NMPC(X(k,:), MPC);  %优化  
            end
        end
      end      
     %F = collectSampleFnc(Xfk, MPC,MPC.w_old,MPC.b_old);
     % 数据预处理,将训练集归一化到[0,1]区间
     % mapminmax为MATLAB自带的归一化函数
%      [Xfk_scale,ps] = mapminmax(Xfk',0,1);
%      Xfk_scale = Xfk_scale';
     % SVM网络训练
     positive = 0;
     negative = 0;
     for k = 1:size(F,1)
         if F(k) == 1
             positive = positive + 1;
         else
             negative = negative + 1;
         end
     end
     model = svmtrain(F, X, '-s 0 -t 2 -c 1 -g 0.07 -b 1');
     positive
     negative
%      % SVM网络预测
%      [predict_label,accuracy,dec_values] = svmpredict(F,Xfk,model); 
     %支持向量索引(Support Vectors Index)
     SVs_idx = model.sv_indices;
     %支持向量特征属性和类别属性
     X_SVs = X(SVs_idx,:);% or use: SVs=full(model.SVs);
     F_SVs = F(SVs_idx);
     %求平面w'x+b = 0的法向量w
     alpha_SVs = model.sv_coef;%实际是a_i*y_i
     w = sum(diag(alpha_SVs)*X_SVs)';%即西瓜书公式(6.9)
     b = -model.rho;
     f_x = X_SVs * w + b;
     plot3(X(F==1,1), X(F==1,2),X(F==1,3), 'b.',...
            X(F==-1,1), X(F==-1,2),X(F==-1,3), 'r.', 'markersize', 14);
%    axis([-1 1 -1 1]);
     hold on;
     xlabel('x_1', 'fontsize', 20);
     ylabel('x_2', 'fontsize', 20);
     h = legend('Feasible', 'Infeasible');
     set(gca, 'fontsize', 20);
     set(h, 'fontsize', 20, 'orientation', 'vertical');
     hold off;
     %判断是否满足循环结束条件
     Stop_sum = 0;
     PNstop = model.totalSV;
     for k = 1:PNstop
         Stop_sum = Stop_sum + abs(w'*X_SVs(k,:)' + b - (MPC.w_old*X_SVs(k,:)' + MPC.b_old));
     end
     if Stop_sum <= epslion * PNstop
         break
     end
     %下一次循环准备
     Xfk = X(F==1,:);
     Uf = U(F==1,:);
     MPC.w_old = w';
     MPC.b_old = b;
     X_SVs_old = X_SVs;
     xfi_num = xfi_num + 1;
     
end
    
%      %稀疏网格数据测试集
%       h = 0.02; % Mesh grid step size
%      [X1,X2] = meshgrid(-1:h:1, -1:h:1);
%      test_label = ones(101*101,1);
%      [~,score] = predict(model,[X1(:), X2(:)]);
%      scoreGrid = reshape(score(:,1), size(X1,1), size(X2,2));

%      %由求网格的等高线画出分类曲面，只有二维能画出来，另一种方法就是在支持向量机的...
%      %核函数空间（三维）中的曲面，将这些三维曲面的网格点映射到二维平面得到点列，既是分类曲线...
%      %三维同理，再高维就画不出来了。而我们这里用简单的等高线来表示，比较方便直观
%      [predict_label,accuracy,dec_values] = svmpredict(dd,[X1(:), X2(:)],model,'-b 1'); %对网格点进行预测，求出网格点的概率值dec_values
%      dec_valuesGrid = reshape(dec_values(:,1), size(X1,1), size(X2,2));%
%      predict_labelGrid = reshape(predict_label(:,1), size(X1,1), size(X2,2));

%% 存储数据
Xf = Xfk;
save('samples_for_ENMPC.mat', 'Xf', 'w', 'b', 'Uf', 'model');
