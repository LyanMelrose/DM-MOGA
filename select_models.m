function [TenTimesFLE,C_index] = select_models(cluster,Network)
global sim_Matrix

C_F=tabulate(cluster);
[~,CF_order_idx]=sort(C_F(:,2),'descend');
TenTimesFLE=[];
for t=1:10
    save_FLE=[];
    save_idx=[];
    for i=1:length(C_F(:,2))
        NowC_idx=find(cluster==CF_order_idx(i));
        save_idx=[save_idx;NowC_idx'];
        fle=FLE(cluster(save_idx)',sim_Matrix(save_idx,save_idx),Network(save_idx,save_idx));
        save_FLE=[save_FLE fle];
    end
    TenTimesFLE=[TenTimesFLE;save_FLE];
end
Avg_TTFLE=mean(TenTimesFLE);
[~,C_index]=max(Avg_TTFLE);
end

