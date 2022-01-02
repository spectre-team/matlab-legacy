function [ thresholds ] = fetch_thresholds( vals, max_components_number )
%FETCH_THRESHOLDS Decomposes parameters into Gaussian peaks and finds
%crossings
%   thresholds = FETCH_THRESHOLDS(vals, max_components_number) - returns
%   cell containing all thresholds. If no crossings are present,
%   no threshold will be returned.

	n = size(vals(:), 1);

	thresholds = [];
    
    max_components_number = double(max_components_number);

    a = cell(1,max_components_number);
    mu = cell(1,max_components_number);
    sig = cell(1,max_components_number);
    TIC = cell(1,max_components_number);
    BIC = Inf*ones(1,max_components_number);
    param = vals(:);
    parfor compnum=1:max_components_number
        [a{compnum},mu{compnum},sig{compnum},TIC{compnum},l_lik]=gaussian_mixture_simple(param,ones(n,1),compnum);
        BIC(compnum) = -2*l_lik+(3*compnum-1)*log(n);
    end

    if isinf(min(BIC))
        return;
    end

    best_compnum_ = find(BIC == min(BIC),1);
    best_compnum = best_compnum_(1);
    best_mu = mu{best_compnum};
    best_sig = sig{best_compnum};
    best_a = a{best_compnum};
    %FIX
    TIC = TIC{best_compnum};

    f = gmm_uborder_fun(best_mu, best_sig, best_a);

    %if there are few components
    if best_compnum>1
        %they are picked pairwise
        for j=1:(best_compnum-1)
            for k=(j+1):best_compnum
                %their upper bound
                g = gmm_uborder_fun(best_mu([j,k]), best_sig([j,k]), best_a([j,k]));
                %the disjunction condition & weight condition is checked
                if g(best_mu(j))==best_a(j)*normpdf(best_mu(j),best_mu(j),best_sig(j)) && ...
                        g(best_mu(k))==best_a(k)*normpdf(best_mu(k),best_mu(k),best_sig(k))
                    %their crossing is found
                    crossing = fminbnd(g,min(best_mu([j,k])),max((best_mu([j,k]))));
                    %and it is checked to be on the upper border
                    if f(crossing)==g(crossing)
                        %if it truly is, then it is saved
                        thresholds = [thresholds; crossing];
                    end
                end
            end
        end
        thresholds = sort(thresholds,'ascend');
    end

end
