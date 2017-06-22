function fll = PSD_fisher_obs_MLE(y, S, gradS, hessS)
% PSD_FISHER (observed) fisher information for Whittle MLE likelihood

d = size(gradS,2);
gradS2 = OuterProduct(gradS);
gradS2 = gradS2./repmat(S.^2, 1,d,d);
fll = hessS./repmat(S, 1,d,d);
fll = gradS2 - repmat((1-y./S), 1,d,d) .* (2*gradS2 - fll);
fll = squeeze(sum(fll));
end

