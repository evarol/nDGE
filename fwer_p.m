function p_out=fwer_p(p_val,alpha)

if nargin<2
    alpha=0.05;
end

[p_vec,idx]=sort(abs(p_val(:)),'ascend','MissingPlacement','last');

for i=1:size(p_vec,1)
    q_thresh=alpha/(length(p_vec)-i+1);
end

p_out(idx)=double(p_vec<=q_thresh);

p_out=reshape(p_out,size(p_val));
