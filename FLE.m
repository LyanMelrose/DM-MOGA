function FLE = FLE(cluster,Module_ori_idx,sim_Matrix,Network)

cnum=length(Module_ori_idx);
FLE=[];
rand_sim=triu(sim_Matrix,1);
rand_sim=reshape(rand_sim,[numel(rand_sim),1]);
rand_sim(rand_sim==0)=[];
rand_sim=rand_sim(randperm(numel(rand_sim),1));
for i=1:cnum
    InC=[];
    CC=find(cluster==Module_ori_idx(i));
    for j=1:length(CC)
        InC=[InC find(Network(CC(j),:)~=0&(cluster==Module_ori_idx(i))')];
    end
    InC=unique(InC);
    InC_sim=triu(sim_Matrix(InC,InC),1);
    InC_sim=reshape(InC_sim,[numel(InC_sim),1]);
    InC_sim(InC_sim==0)=[];
    if isempty(InC_sim)
    else
%         random_sim=InC_sim(randperm(numel(InC_sim),1));
        FLE=[FLE (mean(InC_sim)-rand_sim)];
    end
end

end