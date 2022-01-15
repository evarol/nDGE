function dir_srand=sym_generate_constrained_srand(s1,c1,ntry)
% Syntax:
% dir_srand=dir_generate_srand(s1)
% s1 - the adjacency matrix of a directed network
% ntry - (optional) the number of rewiring steps. If none is given ntry=4*(# of edges in the network)
% Output: dir_srand - the adjacency matrix of a randomized network with the same set of in- and out-degrees as the original one

sd=diag(s1);s1=s1-diag(sd);
cd=diag(c1);c1=c1-diag(cd);
dir_srand=s1;
nrew=0;
[i_srand,j_srand]=find(dir_srand);
[Ne, aux]=size(i_srand);
if (nargin < 3) ntry=4*Ne; end;
for i=1:ntry;
    e2=[];
    while isempty(e2);
        [Ne, aux]=size(i_srand);
        e1=1+floor(Ne*rand);
        v1=i_srand(e1);
        v2=j_srand(e1);
        
        feas_set_1 = find(c1(v1,:));
        feas_set_2 = find(c1(:,v2));
        try;e2=randsample((1:size(i_srand,1))', 1, true,  and(ismember(i_srand,feas_set_2),ismember(j_srand,feas_set_1)));end
    end
    
    v3=i_srand(e2);
    v4=j_srand(e2);
    
    if (v1~=v3)&(v1~=v4)&(v2~=v4)&(v2~=v3);
        if (dir_srand(v1,v4)==0)&(dir_srand(v3,v2)==0);
            
            tmp=dir_srand;
            tmp(v1,v4)=dir_srand(v1,v2);
            tmp(v4,v1)=dir_srand(v2,v1);
            
            tmp(v2,v3)=dir_srand(v3,v4);
            tmp(v3,v2)=dir_srand(v4,v3);
            
            tmp(v1,v2)=0;
            tmp(v2,v1)=0;
            
            tmp(v3,v4)=0;
            tmp(v4,v3)=0;
            
            dir_srand=tmp;clear tmp;
            [i_srand,j_srand]=find(dir_srand);
        end;
    end;
end;


dir_srand=dir_srand+diag(sd);