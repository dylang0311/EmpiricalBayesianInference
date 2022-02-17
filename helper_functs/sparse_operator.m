function L = sparse_operator(params)
% Inputs:
%   - params.N1 = dimension
%   - params.PAORDER = 0 if signal sparse, 1 if edges sparse
%
% Output:
%   - L = sparsifying transform
N1 = params.N1;
N2 = params.N2;

if params.DIM == 1
    switch params.PAORDER
        case {0}
            L = sparse(eye(N1));
        case {1}
            L = sparse(circshift(eye(N1),-1,2) - eye(N1));
    end
elseif params.DIM == 2
    switch params.PAORDER
        case {0}
            L = sparse(eye(N1 * N2));
        case {1}
            T1 = circshift(eye(N1),-1,2) - 2*eye(N1);
            T2 = eye(N1);
            L = sparse(N1*N2,N1*N2);
            L(1:N1,1:N1) = T1;
            L(1:N1,end-N1+1:end) = T2;
            for ii = 1:N2-1
                L(ii*N1+1:(ii+1)*N1,ii*N1+1:(ii+1)*N1) = T1;
                L(ii*N1+1:(ii+1)*N1,(ii-1)*N1+1:ii*N1) = T2;
            end
            L = sparse(L);
    end     
end