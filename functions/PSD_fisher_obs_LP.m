function [fll] = PSD_fisher_obs_LP(S, gradS, hessS, z, MB)
% PSD_FISHER (observed) fisher information using approx log-likelihood
S = reshape(S, [max(size(S)) 1]);
z = reshape(z, [max(size(z)) 1]);

d = size(gradS,2);
gradS2 = OuterProduct(gradS);
gradS2 = gradS2./repmat(S, 1,d,d);
fll = gradS2./repmat(S, 1,d,d) - repmat(((log(S)-z)./S), 1, d, d) .* (gradS2 - hessS);
fll = repmat(MB, length(S), d, d).*fll;
fll = squeeze(sum(fll));
end


