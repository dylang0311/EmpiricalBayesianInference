function L = sparse_operator(params)
% Inputs:
%   - params.N1 = dimension
%   - params.PAORDER = 0 if signal sparse, 1 if edges sparse
%
% Output:
%   - L = sparsifying transform

switch params.PAORDER
    case {0}
        L = eye(params.N1);
    case {1}
        L = circshift(eye(params.N1),-1,2) - eye(params.N1);
end
end