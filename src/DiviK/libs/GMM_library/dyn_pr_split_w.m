%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        main dynamic programming algorithm
%

function [Q,opt_part]=dyn_pr_split_w(data,ygreki,K_gr,aux_mx)

% initialize
Q=zeros(1,K_gr);
N=length(data);
p_opt_idx=zeros(1,N);
p_aux=zeros(1,N);
opt_pals=ones(K_gr,N);
for kk=1:N;
    p_opt_idx(kk)=my_qu_ix_w(data(kk:N),ygreki(kk:N));
end

% aux_mx - already computed

% iterate
for kster=1:K_gr
   for kk=1:N-kster
       for jj=kk+1:N-kster+1
           p_aux(jj)= aux_mx(kk,jj)+p_opt_idx(jj);
       end
       [mm,ix]=min(p_aux(kk+1:N-kster+1));
       p_opt_idx(kk)=mm;
       opt_pals(kster,kk)=kk+ix(1);
   end
   Q(kster)=p_opt_idx(1);
end


% restore optimal decisions
opt_part=zeros(1,K_gr);
opt_part(1)=opt_pals(K_gr,1);
for kster=K_gr-1:-1:1
% 	if K_gr-kster+1<=0
% 		fprintf('K_gr-kster+1<=0\n');
% 	elseif opt_part(K_gr-kster)<=0
% 		fprintf('opt_part(K_gr-kster)<=0\n');
% 	end
	opt_part(K_gr-kster+1)=opt_pals(kster,opt_part(K_gr-kster));
end
