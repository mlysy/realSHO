function [fll] = PSD_hessVec_NLS(S, gradS, hessS, ybinned)
% Calculates the shessian for Normal(S,1) which is "incorrect model" of 
% mean(yPSD) with bin size B
S = reshape(S, [max(size(S)) 1]);
ybinned = reshape(ybinned, [max(size(ybinned)) 1]);
d = size(gradS,2);

Sr = repmat(S, 1, d, d);
ybinnedr = repmat(ybinned, 1, d, d);
fll = ybinnedr.*hessS - (OuterProduct(gradS) + hessS.*Sr);
fll = squeeze(sum(fll));
end

