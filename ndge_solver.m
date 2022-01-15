function [beta_backward,beta_backward_perm,alpha_backward,pval]=ndge_solver(T,C,A0,tf_names,numperm,directed,cutoff_style)

set(groot, 'defaultTextInterpreter','none');
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
vec=@(x)(x(:));

TT=log1p(kron(T,T));
d=size(T,2);
if isempty(A0)
    A=reshape(pinv(C(:)),size(C))-reshape(pinv(1-C(:)),size(C)); % actual contrast matrix -- synapses vs. no synapses
else
    A=reshape(pinv(C(:)),size(C))-reshape(pinv(A0(:)-C(:)),size(C)); % actual contrast matrix -- synapses vs. adjacency
end

beta_backward = reshape(vec(A)'*TT,[size(T,2) size(T,2)]);
beta_backward_perm(:,:,1)=beta_backward;
globalTic=tic;
for t=2:numperm
    if isempty(A0)
        if directed==1
            Cperm=dir_generate_srand(C);
        else
            Cperm=sym_generate_srand(C);
        end
        A=reshape(pinv(Cperm(:)),size(C))-reshape(pinv(1-Cperm(:)),size(C));
    else
        if directed==1
            Cperm=dir_generate_constrained_srand(C,A0);
        else
            Cperm=sym_generate_constrained_srand(C,A0);
        end
        A=reshape(pinv(Cperm(:)),size(C))-reshape(pinv(A0(:)-Cperm(:)),size(C));
    end
    subplot(1,2,1)
    imagesc(Cperm);drawnow
    subplot(1,2,2)
    imagesc(A);drawnow
    beta_backward_perm(:,:,t) = reshape(vec(A)'*TT,[size(T,2) size(T,2)]);
    clc;disp(['Permutation ' num2str(t) '/' num2str(numperm) '. Time left: ' num2str((toc(globalTic)/t)*(numperm-t)) ' seconds.']);
end
close all



alpha_backward = (beta_backward - mean(beta_backward_perm,3))./std(beta_backward_perm,[],3);
pval = 2-2*normcdf(abs(alpha_backward),0,1);
alpha_backward(isnan(alpha_backward))=0;

if strcmpi(cutoff_style,'bonferroni')
    cutoff=icdf('normal',1-0.05/d,0,1);
elseif strcmpi(cutoff_style,'fwer')
    fwer_cutoff=fwer_p(pval,0.05);
    cutoff=min(abs(alpha_backward(fwer_cutoff==1)));
end

%% plot 1
figure('units','normalized','outerposition',[0 0 1 1])
imagesc((alpha_backward.*(abs(alpha_backward)>=cutoff)),[-max(abs(alpha_backward(:))) max(abs(alpha_backward(:)))]);colormap(othercolor('BuDRd_12'));colorbar;set(gca,'Xtick',(1:size(T,2)),'Xticklabel',tf_names,'xticklabelrotation',90,'Ytick',(1:size(T,2)),'Yticklabel',tf_names);
set(gca,'FontWeight','bold','TickLength',[0 0]);
idx=find(abs(diag(alpha_backward))>=cutoff);
for i=1:length(idx)
    text(idx(i)+1,idx(i),tf_names{idx(i)},'Color','k','FontWeight','bold');
end
ylabel('Pre-synaptic partner gene');
xlabel('Post-synaptic partner gene');
title('Heterotypic interactions, Significant Hits');
%% plot 2
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(eye(d).*(alpha_backward.*(abs(alpha_backward)>=cutoff)),[-max(abs(alpha_backward(:))) max(abs(alpha_backward(:)))]);colormap(othercolor('BuDRd_12'));colorbar;set(gca,'Xtick',(1:size(T,2)),'Xticklabel',tf_names,'xticklabelrotation',90,'Ytick',(1:size(T,2)),'Yticklabel',tf_names);
set(gca,'FontWeight','bold','TickLength',[0 0]);
idx=find(abs(diag(alpha_backward))>=cutoff);
for i=1:length(idx)
    text(idx(i)+1,idx(i),tf_names{idx(i)},'Color','k','FontWeight','bold');
end
ylabel('Pre-synaptic partner gene');
xlabel('Post-synaptic partner gene');
title('Homomeric interactions, Significant Hits');

%% plot 3
figure('units','normalized','outerposition',[0 0 1 1])
beta_demeaned=beta_backward-mean(beta_backward_perm,3);
plot(vec(diag(beta_backward-mean(beta_backward_perm,3))),vec(diag(abs(alpha_backward))),'ko','MarkerFaceColor','k','MarkerSize',3);
hold on
for i=1:length(idx)
    if beta_demeaned(idx(i),idx(i))>0
        plot(beta_demeaned(idx(i),idx(i)),abs(alpha_backward(idx(i),idx(i))),'ko','MarkerFaceColor','r','MarkerSize',10);
        text(beta_demeaned(idx(i),idx(i))+max(abs(beta_demeaned(:)))/100,abs(alpha_backward(idx(i),idx(i))),tf_names{idx(i)},'Color','k','FontWeight','bold','FontSize',14);
    else
        plot(beta_demeaned(idx(i),idx(i)),abs(alpha_backward(idx(i),idx(i))),'ko','MarkerFaceColor','b','MarkerSize',10);
        text(beta_demeaned(idx(i),idx(i))+max(abs(beta_demeaned(:)))/100,abs(alpha_backward(idx(i),idx(i))),tf_names{idx(i)},'Color','k','FontWeight','bold','FontSize',14);
    end
end

xlabel('log-fold-change','FontSize',14);
ylabel('log-p-value','FontSize',14);
grid on
set(gca,'FontWeight','bold');
title('Homomeric interactions, Significant Hits','FontSize',14);

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor','r','MarkerSize',8);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor','b','MarkerSize',8);
h(3) = plot(NaN,NaN,'ko','MarkerFaceColor','k','MarkerSize',3);
legend(h, {'Significant synaptic association','Significant synaptic inhibition','Not significant'},'FontSize',14);set(gcf,'Color','w');
%% plot 4
figure('units','normalized','outerposition',[0 0 1 1])
plot(vec(beta_backward-mean(beta_backward_perm,3)),vec(abs(alpha_backward)),'ko','MarkerFaceColor','k','MarkerSize',3);
hold on
sig_pos=find(alpha_backward>=cutoff);
sig_neg=find(alpha_backward<=-cutoff);
plot(beta_demeaned(sig_pos),abs(alpha_backward(sig_pos)),'ko','MarkerFaceColor','r','MarkerSize',10);
plot(beta_demeaned(sig_neg),abs(alpha_backward(sig_neg)),'ko','MarkerFaceColor','b','MarkerSize',10);
[I,J]=find(abs(alpha_backward)>=cutoff);
for i=1:size(I,1)
    %     if I(i)>=J(i)
    text(beta_demeaned(I(i),J(i))+max(abs(beta_demeaned(:)))/100,abs(alpha_backward(I(i),J(i))),[tf_names{I(i)} ' x ' tf_names{J(i)}],'Color','k','FontWeight','bold','FontSize',14);
    %     end
end
xlabel('log-fold-change','FontSize',14);
ylabel('log-p-value','FontSize',14);
title('Heterotypic interactions, Significant Hits','FontSize',14);
grid on
set(gca,'FontWeight','bold');
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor','r','MarkerSize',8);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor','b','MarkerSize',8);
h(3) = plot(NaN,NaN,'ko','MarkerFaceColor','k','MarkerSize',3);
legend(h, {'Significant synaptic association','Significant synaptic inhibition','Not significant'},'FontSize',14);set(gcf,'Color','w');