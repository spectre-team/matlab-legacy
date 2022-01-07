function [ vals ] = normpdfs(x, mu, sig, amp)
n = length(mu);
vals = NaN(n, length(x));
for i=1:n
    vals(i, :) = amp(i) * normpdf(x, mu(i), sig(i));
end
end
