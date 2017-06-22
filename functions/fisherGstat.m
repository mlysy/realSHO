function prob = fisherGstat(a, q, logSort)
% FISHERGSTAT Let U = U0 < ... < Uq, where U1,...,U{q-1} are the _order
% statistics_ of q-1 iid Uniforms, and U0 = 0, Uq = 1.  Let Mq = max(diff(U)).
% Then this function calculates P(Mq > a).
% It turns out that if yPSD with no decimation is applied to white noise,
% then Mq = max(yPSD)/sum(yPSD), so you can use this to test the hypothesis
% of white noise, if Mq is too high you reject H0.

if nargin < 3, logSort = true; end

% Maximum number of values to calculate
maxq = min(q, floor(1/a));
maxq2 = ceil(maxq/2);

% Calculate the terms on the log scale
logt = zeros(1,maxq2*2);
ind = 1:maxq;
logt(ind) = gammaln(q+1) - gammaln(q-ind+1) -gammaln(ind+1) + (q-1)*log(1-ind*a);
if maxq < 2*maxq2, logt(end) = -Inf; end

% Group the +ve and -ve terms by two to calculate the sum
logt = reshape(logt, 2, ceil(maxq/2));

% If this is true try to group terms in the sum to minimize numerical instability
if logSort, logt = sort(logt,2); end
maxlt = max(logt);
logt = exp(logt - repmat(maxlt,2,1));
logt = logt(1,:) - logt(2,:);

% If both terms are minuscule might get Inf-Inf
logt(isnan(logt)) = 0;

% Add the remaining terms
sgt = sign(logt);
logt = log(abs(logt)) + maxlt;
mx = max(logt);
prob = sum(sgt .* exp(logt - mx));
prob = exp(log(prob) + mx);
end

