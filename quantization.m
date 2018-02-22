function [Mq,part,codebook,mn,mx]=quantization(M,bit)
% quantization
% Inputs: 
%         M: a n1 x n2 x ... x nN array of the original values
%         bit = 1 bit of quantization
% Output: 
%         Mq: a n1 x n2 x ... x nN array of quantized values
%         part: scalar quantization bin boundary
%         codebook: a 1 x 2 vector of quantization levels
%         mn: minimum scalar of M
%         mx: maximum scalar of M

mn = min(M(:));
mx = max(M(:));

part = mean(M(:));
codebook = 1:(2^bit);

[~,Mv] = quantiz(M(:),part,codebook);
Mq = reshape(Mv,size(M));
end