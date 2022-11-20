% Effect Size
% import GSE19804_Sanger_ENTREZ_ID----Entrez Id are string matrix, others
% are matrix
% import label of samples(_ENid_BP_GSE19804_Final.label),
% 1--cancer,0--normal
label=logical(label);
n1=sum(label);
n2=length(label)-n1;
s1=std(exp_data(:,label),0,2);
s2=std(exp_data(:,~label),0,2);
Sp=sqrt(((n1-1).*(s1.^2)+(n2-1).*(s2.^2))./(n1+n2-2));
g_score=(mean(exp_data(:,label),2)-mean(exp_data(:,~label),2))./Sp;

fpath='DEG_result_Effect_Size_GSE19188.txt';
fid=fopen(fpath,'wt');%写入文件路径
out=[Entrez_Id g_score];                            %input_matrix为待输出矩阵
[m,n]=size(out);
for i=1:m
    for j=1:n
        x=out(i,j);
        if j~=n
            fprintf(fid,'%s\t',x);
        else
            fprintf(fid,'%s\n',x);
        end
    end
end
%撤销文件句柄
fclose(fid);

save DEG_result_Effect_Size_GSE19188.mat g_score Entrez_Id

% % 根据 g_score 的降幂排序来对 Entrez_Id 进行排序
% % 对两个数据集都排序处理，结果加数据集编号前缀来保存
% [sorted_g_score,ori_id]=sort(g_score,'descend');
% sorted_Entrez_Id=Entrez_Id(ori_id,1);

% 两数据集都只选绝对值大于0.8的DEGs做重叠,
% 正值为癌症样本上调，负值为癌症样本下调
GSE19804_DEGs_ENid=GSE19804_ENid(abs(GSE19804_g_score)>0.8,1);
GSE19188_DEGs_ENid=GSE19188_ENid(abs(GSE19188_g_score)>0.8,1);
% 观察两个数据集的 DEGs 的重叠部分
overlap_idx=ismember(GSE19804_DEGs_ENid,GSE19188_DEGs_ENid);
overlap_ENid_Effect_Size=GSE19804_DEGs_ENid(overlap_idx);
save overlap_ENid_Effect_Size.mat overlap_ENid_Effect_Size;
% 观察Effect Size所得重叠基因和limma所得重叠基因间是否重叠
oo_DEGs=ismember(overlap_ENid_Effect_Size,overlap_ENid_limma);
all(oo_DEGs)

% 检查上调和下调基因在模块中比例是否变高了
GSE19804_M1_UPgene=sum(ismember(GSE19804_M1genes,GSE19804_UPgene));
GSE19804_M1_DOWNgene=sum(ismember(GSE19804_M1genes,GSE19804_DOWNgene));
GSE19188_M1_UPgene=sum(ismember(GSE19188_M1genes,GSE19188_UPgene));
GSE19188_M1_DOWNgene=sum(ismember(GSE19188_M1genes,GSE19188_DOWNgene));