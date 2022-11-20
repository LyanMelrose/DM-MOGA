function idx = DBIndex(dist,CC_idx,sampleCI)

%dist=dist,sampleCI=cl,CNum=NCLUST,CCIndex=icl������������ԭ���ݼ��е��±�
%����DB,Iָ��,
% %DMΪԭʼ���ݾ���
% global DM;

%�������
[ci,~]=max(sampleCI);

%�ݴ��ж�
if ci>1&&ci<length(sampleCI)

    for i=1:ci
        %�ҵ����ڴ�i���±�ֵ
        y=find(sampleCI==i);
        %�õ�������ж��ٸ�����
        nx=length(y);
        %ȫ����0�����
        Sumi=sum(dist(CC_idx(i),y));
        %����0
        Si(i)=Sumi/nx;
    end
    
    %��ʼ��Rֵ
    R=zeros(ci,ci);
    
    %ͳ�Ƽ�����ؼ��Rֵ
    for i=1:ci-1
        for j=i+1:ci
            %����i��j���Rֵ
            if dist(CC_idx(i),CC_idx(j))==0
                R(i,j)=0;
            else
                R(i,j)=(Si(i)+Si(j))/dist(CC_idx(i),CC_idx(j));
            end
            R(j,i)=R(i,j);
        end
        %����õ���i��DBָ��ֵ
        DB(i)=max(R(i,:))/ci;
    end
    
    %����õ�������
    for i=1:ci
        for j=i:ci
            dk(i,j)=dist(CC_idx(i),CC_idx(j));
        end
    end
    
    %�õ�dk�����ֵ
    DK=max(max(dk));
    %����õ�I_index��ֵ
%     I_index=((E1*DK)/(ci*EK))^2;
    
    %����Silָ��
%     R=silhouette(DM,sampleCI, 'euclidean');
    %��Silָ���ƽ��ֵ
%     Sil=mean(R);
    
%     idx=Sil;
    %�ܵ�DBָ��
    idx=max(DB);

else
    idx=0;
end

end
