function [ f] = gmm_uborder_fun( mu, sig, amp )
	if sum(size(mu)==size(sig))<length(size(mu)) || sum(size(mu)==size(amp))<length(size(mu))
		error('Inconsistent sizes of input arguments');
	end
	f = @(x) max(normpdfs(x, mu(:), sig(:), amp(:)), [], 1);
	
end