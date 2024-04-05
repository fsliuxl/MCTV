function xrec = mctv_mc(Y,omega,opts)
% Solve min_X sum_{i=1}^2||T_i(X)||_{w,*}, s. t. P_Omega(X) = P_Omega(X_0)
% Input: Y
%% opts
X      = opts.X;% used for calculating the quality indices
mu     = opts.mu;
rho    = opts.rho;
MaxIter= opts.MaxIter;
tol    = opts.tol;
%% Initialization
dim    = size(Y);
X0     = zeros(dim);
E      = X0;
Z1     = E; Z2 = E; DualZ1 = E; DualZ2 = E; DualZ3 = E;
%% Main Loops
for k = 1:MaxIter
    tic;
    % Update X
    rhs = reshape(diff_x(Z1 - DualZ1/mu,dim) + diff_y(Z2 - DualZ2/mu,dim),dim) - E + Y - DualZ3/mu;
    x   = myPCG_mc(X0(:),rhs(:),dim);
    X1  = reshape(x,dim);
    % Update E
    E = Y - X1 - DualZ3/mu; E(omega) = 0;
    % Update Z1
    [U,S,V] = svd(reshape(diff_x(X1,dim),dim) + 1/mu*DualZ1,'econ');
    Z1 = U*MySoftTh(S,1/mu)*(V');
    % Update Z2
    [U,S,V] = svd(reshape(diff_y(X1,dim),dim) + 1/mu*DualZ2,'econ');
    Z2 = U*MySoftTh(S,1/mu)*(V');
    % Update DualZ1 and DualZ2 and DualZ3
    DualZ1 = DualZ1 + mu * (reshape(diff_x(X1,dim),dim) - Z1);
    DualZ2 = DualZ2 + mu * (reshape(diff_y(X1,dim),dim) - Z2);
    DualZ3 = DualZ3 + mu * (X1 + E - Y);
    % Calculate
    err = norm(X0-X1,'fro')/max(norm(X0,'fro'),1);
    err = min(err,1);
    [psnr1,ssim1] = msqia(reshape(X,dim),reshape(X1,dim));
    fprintf('iter = %g, cost time = %.4f, mu=%.4f, RelErr = %.6f, PSNR = %.4f, SSIM = %.4f\n',k, toc, mu,err,psnr1,ssim1);
    if err <= tol
        fprintf('>>>>>>>>>>>>>>>> ADMM convergenced >>>>>>>>>>>>>>>>\n');
        break;
    end
    X0 = X1;
    mu = min(mu * rho, 1e8);
end
xrec = reshape(X1,dim);
end
%% This is soft thresholding operation
function X= MySoftTh(B,lambda)
X=sign(B).*max(0,abs(B)-lambda);
end