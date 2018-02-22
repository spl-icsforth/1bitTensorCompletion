function W = Fold(W, dim, i)
% Fold a matricization mode into a tensor
dim = circshift(dim, [1-i, 1-i]);
W = shiftdim(reshape(W, dim), length(dim)+1-i);
end