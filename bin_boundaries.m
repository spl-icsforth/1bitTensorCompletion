function [U,L] = bin_boundaries(Mq,part,codebook,mn,mx)
% Compute the upper and lower quantization bin boundaries of each measurement
% Input: 
%         Mq: a n1 x n2 x ... x nN array of quantized values
%         part: scalar quantization bin boundary
%         codebook: a 1 x 2 vector of quantization levels
%         mn: minimum scalar of M
%         mx: maximum scalar of M
%
% Output:
%         U: a n1 x n2 x ... x nN array of the upper bin boundaries of Mq
%         L: a n1 x n2 x ... x nN array of the lower bin boundaries of Mq

Nway = size(Mq);
data = Mq(:);
U = zeros(size(data));
L = zeros(size(data));
for j = 1:length(data)
    for k = 1:length(part)
        if (data(j) == codebook(k))
            if (k == 1)
                L(j) = mn;
                U(j) = part(k);
               break;
            else
                L(j) = part(k-1);
                U(j) = part(k);
                break;
            end
        elseif (k == length(part))
            L(j) = part(k);
            U(j) = mx;
            break;
        end
    end
end
U = reshape(U,Nway);
L = reshape(L,Nway);
end