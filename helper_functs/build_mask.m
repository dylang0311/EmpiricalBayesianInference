function [T,zeroColumns] = build_mask(m,params)
L = sparse_operator(params);
T = m(:).*L;
zeroColumns = [];
% 
% T(~any(T,2),:) = [];
% 
% zeroColumns = ~any(T,1);
% T(:,zeroColumns) = [];
end