function F = collectSampleFnc(xk, MPC,w_svm,b_svm)
    Ns = size(xk,1);
    F = zeros(Ns, 1);   % allocating memory for feasibility information
    for k = 1:Ns
        fprintf('\nIteration: %d of %d', k, Ns);
    %     [U(k), F(k)] = find_optimal_NMPC(X(k,:), MPC);
        [F(k)] = find_optimal_NMPC(xk(k,:), MPC,w_svm,b_svm);
    end