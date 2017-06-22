function [B,grad] = PSD_gradVec_NLS(S, gradS, ybinned)
% Calculates the first derivative of the log-likelihood (score) for
% Normal(S,1) which is distribution of mean(yPSD) with bin size B
S = reshape(S, [max(size(S)) 1]);
ybinned = reshape(ybinned, [max(size(ybinned)) 1]);
d = size(gradS,2);

Sr = repmat(S, 1, d);
ybinnedr = repmat(ybinned, 1, d);
g = (ybinnedr - Sr).*gradS;
grad = sum(g,1);
B = g' * g;
end

