function [Lrec, Srec] = mctv_rpca(Y,opts)
% Solve min_X sum_{i=1}^2||T_i(X)||_{w,*} + lambda*|S|, s. t. Y = X + S
% Input: Y
%% opts
X      = opts.X;% used for calculating the quality indices
mu     = opts.mu;
rho    = opts.rho;
MaxIter= opts.MaxIter;
tol    = opts.tol;
lambda = opts.lambda;

%% Initialization
dim    = size(Y);
L0     = Y;
S0      = L0;
Z1     = S0; Z2 = S0; DualZ1 = S0; DualZ2 = S0; DualZ3 = S0;
%% Main Loops
for k = 1:MaxIter
    tic;
    % Update X
    rhs = reshape(diff_x(Z1 - DualZ1/mu,dim) + diff_y(Z2 - DualZ2/mu,dim),dim) - S0 + Y - DualZ3/mu;
    x   = myPCG_rpca(L0(:),rhs(:),dim);
    L1  = reshape(x,dim);
    S1 = MySoftTh(Y-L1-DualZ3/mu, lambda/mu);
    % Update Z1
    [U,S,V] = svd(reshape(diff_x(L1,dim),dim) + 1/mu*DualZ1,'econ');
    Z1 = U*MySoftTh(S,1/mu)*(V');
    % Update Z2
    [U,S,V] = svd(reshape(diff_y(L1,dim),dim) + 1/mu*DualZ2,'econ');
    Z2 = U*MySoftTh(S,1/mu)*(V');
    % Update DualZ1 and DualZ2 and DualZ3
    DualZ1 = DualZ1 + mu * (reshape(diff_x(L1,dim),dim) - Z1);
    DualZ2 = DualZ2 + mu * (reshape(diff_y(L1,dim),dim) - Z2);
    DualZ3 = DualZ3 + mu * (L1 + S1 - Y);
    % Calculate
    err_L = norm(L0-L1,'fro')/max(norm(L0,'fro'),1);
    err_L = min(err_L,1);
    err_S = norm(S0-S1,'fro')/max(norm(S0,'fro'),1);
    err_S = min(err_S,1);
    [psnr1,ssim1] = msqia(reshape(X,dim),reshape(L1,dim));
    fprintf(['iter = %g, cost time = %.4f, Rel_L = %.6f, Rel_S = %.6f,'...
        'PSNR = %.4f, SSIM = %.4f\n'],k, toc,err_L,err_S,psnr1,ssim1);
    if err_L <= tol
        fprintf('>>>>>>>>>>>>>>>> ADMM convergenced >>>>>>>>>>>>>>>>\n');
        break;
    end
    L0 = L1;
    S0 = S1;
    mu = min(mu * rho, 1e8);
end
Lrec = reshape(L1,dim);
Srec = reshape(S1,dim);
end
%% This is soft thresholding operation
function X= MySoftTh(B,lambda)
X=sign(B).*max(0,abs(B)-lambda);
end
