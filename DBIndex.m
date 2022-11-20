function idx = DBIndex(dist,CC_idx,sampleCI)

%dist=dist,sampleCI=cl,CNum=NCLUST,CCIndex=icl即聚类中心在原数据集中的下标
%包含DB,I指数,
% %DM为原始数据矩阵
% global DM;

%聚类个数
[ci,~]=max(sampleCI);

%容错判断
if ci>1&&ci<length(sampleCI)

    for i=1:ci
        %找到属于簇i的下标值
        y=find(sampleCI==i);
        %得到这个簇有多少个样本
        nx=length(y);
        %全等于0的情况
        Sumi=sum(dist(CC_idx(i),y));
        %等于0
        Si(i)=Sumi/nx;
    end
    
    %初始化R值
    R=zeros(ci,ci);
    
    %统计计算各簇间的R值
    for i=1:ci-1
        for j=i+1:ci
            %计算i和j间的R值
            if dist(CC_idx(i),CC_idx(j))==0
                R(i,j)=0;
            else
                R(i,j)=(Si(i)+Si(j))/dist(CC_idx(i),CC_idx(j));
            end
            R(j,i)=R(i,j);
        end
        %计算得到簇i的DB指数值
        DB(i)=max(R(i,:))/ci;
    end
    
    %计算得到类间距离
    for i=1:ci
        for j=i:ci
            dk(i,j)=dist(CC_idx(i),CC_idx(j));
        end
    end
    
    %得到dk的最大值
    DK=max(max(dk));
    %计算得到I_index的值
%     I_index=((E1*DK)/(ci*EK))^2;
    
    %计算Sil指标
%     R=silhouette(DM,sampleCI, 'euclidean');
    %求Sil指标的平均值
%     Sil=mean(R);
    
%     idx=Sil;
    %总的DB指数
    idx=max(DB);

else
    idx=0;
end

end
