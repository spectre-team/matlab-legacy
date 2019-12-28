% EM iterations 
function [ppsort,musort,sigsort,l_lik] = g_mix_est_fast_lik(raw_sample,KS,muv,sigv,pp)
% input data:
% sample - vector of observations
% muv,sigv,pp - initial values of mixture parameters
% output:
% musort,sigsort,ppsort - estimated mixture parameters
% l_lik - log likeliood

min_sg=0.1;

raw_sigvoc=max(sigv.^2, min_sg^2);

sample = raw_sample(:)';
mivoc = muv(:);
sigvoc = raw_sigvoc(:);
ppoc = pp(:);



N=length(sample);

% OK oceniamy iteracyjnie wg wzorow z artykulu Bilmsa
pssmac=zeros(KS,N);
change=1;
rescue_iterations_count = 0;
unreachable_number_of_iterations = 10000;

while change > 1.5e-4 && rescue_iterations_count < unreachable_number_of_iterations
   rescue_iterations_count = rescue_iterations_count + 1;
   oldppoc=ppoc;
   oldsigvoc=sigvoc;
   
   % lower limit for component weights
   ppoc=max(ppoc,0.001);
     
   for kskla=1:KS
      pssmac(kskla,:)=ppoc(kskla)*normpdf(sample,mivoc(kskla),sqrt(sigvoc(kskla)));
   end
   psummac=ones(KS,1)*sum(pssmac,1);
   pskmac=pssmac./psummac;
   ppp=sum(pskmac,2);
   ppoc=ppp/N;
   mivoc=pskmac*sample';
   mivoc=mivoc./ppp;
   sigmac=(ones(KS,1)*sample-mivoc*ones(1,N)).*((ones(KS,1)*sample-mivoc*ones(1,N)));
   for kkk=1:KS
      % lower limit for component variances 
      sigvoc(kkk)=max([pskmac(kkk,:)*sigmac(kkk,:)' min_sg^2]);
   end
   sigvoc=sigvoc./ppp;
 
   %
   change=sum(abs(ppoc-oldppoc))+sum(abs(sigvoc-oldsigvoc))/KS;
end

% compute likelihood
l_lik=sum(log(sum(pssmac,1))); 

% sort estimates
[musort,isort]=sort(mivoc);
sigsort=sqrt(sigvoc(isort));
ppsort=ppoc(isort);



