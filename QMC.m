function Mhat = QMC(Y,idx,U,L,model)
% Quantized Matrix Completion
% Inputs: 
%         Y: d1 x d2 measurement matrix
%         idx: p x 1 vector of indices of the observations
%         U: d1 x d2 matrix of the upper bin boundaries of each measurement
%         L: d1 x d2 matrix of the lower bin boundaries of each measurement
%         model: 1 for logistic, 2 for probit
% Output: 
%         Mhat the d1 x d2 estimated matrix


%  parameters
itter = 50;
threshold = 1e-7;


[d1, d2] = size(Y);
Z1 = Y;
Z2 = zeros(size(Y));
df = zeros(size(Y));
    
if model == 1
    L_log = 1.0/4.0; % Lipschitz
    s = 1.0 / L_log; % step size

    i = 1;
    error = 0;
    while ((i <= itter)&&(error >= threshold)) || (i == 1)
        if i~=1
            Z1 = Z2;
        end
        for k = 1:length(idx)
            if Z1(idx(k)) > U(idx(k))
                Z1(idx(k)) = U(idx(k));
            end
            if Z1(idx(k)) < L(idx(k))
                Z1(idx(k)) = L(idx(k));
            end
        end
        fprintf('Iteration: %d\n',i);

        % Reduce the objective function f(Z)
        df(idx) = (f_log_prime(L(idx)-Z1(idx))-f_log_prime(U(idx)-Z1(idx)))./(f_log(U(idx)-Z1(idx))-f_log(L(idx)-Z1(idx)));

        Zhat = Z1(:)-s.*df(:);
        Zhat = reshape(Zhat,d1,d2);
        
        % Impose low-rankness on Zhat
        %sparse
		[I,J] = ind2sub([d1,d2], 1:d1*d2);
		V = Zhat(:);
		D = spconvert([I',J',V; d1,d2,0]);
        %matrix completion algorithm - ALM method
		[B,~,~] = inexact_alm_mc(D, 1e-7,17);
		Z2 = B.U*B.V';
        
        Z2 = round(Z2);
        
        i = i+1;
        error = norm(Z2-Z1,'fro')/norm(Z2,'fro');
    end
else
    L_pro = 1.0; % Lipschitz
    s = 1.0 / L_pro; % step size
    % mean value and standart deviation for the inverse probit link function
    mu = 0;
    sigma = 1;

    i = 1;
    error = 0;
    while ((i <= itter)&&(error >= threshold)) || (i == 1)
        if i ~= 1
            Z1 = Z2;
        end
        for k = 1:length(idx)
            if Z1(idx(k)) > U(idx(k))
                Z1(idx(k)) = U(idx(k));
            end
            if Z1(idx(k)) < L(idx(k))
                Z1(idx(k)) = L(idx(k));
            end
        end
        fprintf('Iteration: %d\n',i);
        
        % Reduce the objective function f(Z)
        df(idx) = (f_pro_prime(L(idx)-Z1(idx),mu,sigma)-f_pro_prime(U(idx)-Z1(idx),mu,sigma))./(f_pro(U(idx)-Z1(idx),mu,sigma)-f_pro(L(idx)-Z1(idx),mu,sigma));

        Zhat = Z1(:)-s.*df(:);
        Zhat = reshape(Zhat,d1,d2);
        
        % Impose low-rankness on Zhat
        %sparse
		[I,J] = ind2sub([d1,d2], 1:d1*d2);
		V = Zhat(:);
		D = spconvert([I',J',V; d1,d2,0]);
        %matrix completion algorithm - ALM method
		[B,~,~] = inexact_alm_mc(D, 1e-7,17);
		Z2 = B.U*B.V';
        
        Z2 = round(Z2);

        i = i+1;
        error = norm(Z2-Z1,'fro')/norm(Z2,'fro');
    end
end
Mhat = Z2;
end