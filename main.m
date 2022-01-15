clear all
clc
close all

%% Inputs from user

%INSERT GENES THAT SHOULD BE INCLUDED
wanted_gene_names={'unc-42','lad-2','ncam-1','nlg-1','oig-1','rig-1','rig-3','rig-5','rig-6','syg-1','unc-6','zig-5','zig-8'};
wanted_neurons=[];%{'ASH','AIB','AVH','AVK','AVD','SAA','AVA','RIV','RMF','AVE','AVB','RMH','SMD','SB','RMD_DV','RMD_LR'}; %%e.g. UNC-42 subcircuit - leave empty as [] to include all
numperm=1000; % number of permutations to run nDGE
binary_connectome=1; %switch to 0 or 1 if binary connectome analysis should be run
binary_gene_expr=1; %switch to 0 or 1 if binary gene expr. analysis should be run
directed_analysis=1; %treat the connectome as directed or symmetric?
cutoff_style='fwer'; %FDR cutoff - options: 'bonferroni', 'fwer'
%% For NDGE we need 3 things, gene expression, connectivity, and adjacency (if available)

%% importing gene expression
disp('Choose gene expression file (genes x neurons)');
[filename, pathname]=uigetfile('*','Choose gene expression file (genes x neurons)');
gene_expr_table=readtable([pathname filename]);
gene_expr_neurons=gene_expr_table.Properties.VariableNames(2:end);
gene_expr_genes=table2array(gene_expr_table(:,1));
gene_expr_matrix=table2array(gene_expr_table(:,2:end));

%% importing connectome
disp('Choose connectivity file');
[filename, pathname]=uigetfile('*','Choose connectivity file');

connectivity_table=readtable([pathname filename]);
connectivity_matrix = table2array(connectivity_table);

%% importing connectome row/column names
disp('Choose connectivity neuron names file');
[filename, pathname]=uigetfile('*','Choose connectivity neuron names file');

connectivity_names_table=readtable([pathname filename],'readvariablenames',false);
connectivity_neuron_names = table2array(connectivity_names_table);


try
    %% importing adjacency
    disp('Choose adjacency file');
    [filename, pathname]=uigetfile('*','Choose adjacency file');
    adjacency_table=readtable([pathname filename]);
    adjacency_matrix = table2array(adjacency_table);
    %% importing adjacency row/column names
    disp('Choose connectivity neuron names file');
    [filename, pathname]=uigetfile('*','Choose connectivity neuron names file');
    adjacency_names_table=readtable([pathname filename],'readvariablenames',false);
    adjacency_neuron_names = table2array(adjacency_names_table);
    
catch
    warning('no adjacency provided - will proceed with only connectome');
end

%% matching common neurons across genomes and connectomes and wanted neurons
if ~exist('adjacency_matrix')
    if ~isempty(wanted_neurons)
        common_neurons=intersect(intersect(gene_expr_neurons,connectivity_neuron_names),wanted_neurons);
    else
        common_neurons=intersect(gene_expr_neurons,connectivity_neuron_names);
    end
    [a,b]=ismember(common_neurons,connectivity_neuron_names);
    C=connectivity_matrix(b(a==1),b(a==1));
    [a,b]=ismember(common_neurons,gene_expr_neurons);
    T=gene_expr_matrix(:,b(a==1));
else
    if ~isempty(wanted_neurons)
        common_neurons=intersect(intersect(intersect(gene_expr_neurons,connectivity_neuron_names),adjacency_neuron_names),wanted_neurons);
    else
        common_neurons=intersect(intersect(gene_expr_neurons,connectivity_neuron_names),adjacency_neuron_names);
    end
    [a,b]=ismember(common_neurons,connectivity_neuron_names);
    C=connectivity_matrix(b(a==1),b(a==1));
    [a,b]=ismember(common_neurons,adjacency_neuron_names);
    A=adjacency_matrix(b(a==1),b(a==1));
    [a,b]=ismember(common_neurons,gene_expr_neurons);
    T=gene_expr_matrix(:,b(a==1));
end


%% subselecting wanted genes
[a,b]=ismember(wanted_gene_names,gene_expr_genes);
T=T(b(a==1),:);
gene_names=gene_expr_genes(b(a==1));


%% binarizing connectomes and genomes if necessary
if binary_connectome==1
    C_input=double(C>0);
    A_input=[];
    if exist('A')
        A_input=double(A>0);
    end
else
    C_input=C;
    A_input=[];
    if exist('A')
        A_input=A;
    end
end


if binary_gene_expr==1
    T_input=double(T>0);
else
    T_input=T;
end

%% main nDGE routine
[beta_backward,beta_backward_perm,alpha_backward,pval]=ndge_solver(T_input',C_input,A_input,gene_names,numperm,directed_analysis,cutoff_style);




%% writing outputs to excel
effect_size_table = array2table((beta_backward - mean(beta_backward_perm,3)));
effect_size_table.Properties.RowNames=matlab.lang.makeValidName(gene_names);
effect_size_table.Properties.VariableNames=matlab.lang.makeValidName(gene_names);

pval_table = array2table(2-2*normcdf(abs(alpha_backward),0,1));
pval_table.Properties.RowNames=matlab.lang.makeValidName(gene_names);
pval_table.Properties.VariableNames=matlab.lang.makeValidName(gene_names);


mkdir('output');
writetable(effect_size_table,'./output/effect_sizes.xlsx','WriteRowNames',true,'WriteVariableNames',true);
writetable(pval_table,'./output/pvals.xlsx','WriteRowNames',true,'WriteVariableNames',true);

