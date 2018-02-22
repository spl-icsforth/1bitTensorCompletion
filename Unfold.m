function W = Unfold(W, dim, i)
% Unfold the tensor into a matricization mode
W = reshape(shiftdim(W,i-1), dim(i), []);

end
