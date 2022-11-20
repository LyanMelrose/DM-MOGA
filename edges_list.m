function [et] = edges_list(A,V,Node3)
%EDGELIST Generate the edges list for each node.

for i=1:V
    for j=1:V
        if i~=j
            if A(i,i)==0
                m=1;
            else
                m=single(A(i,j))/single(A(i,i));    %局部社团 i 的外部边权与内部边权的比值。
            end
            if A(j,j)==0
                k=1;
            else
                k=single(A(i,j))/single(A(j,j));    %局部社团 j 的外部边权与内部边权的比值。
            end
            if m<k
                m=k;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if m<1/4
%                 A(i,j)=0;
%             end
        end
    end
end

%%
et = [];
for i = 1:V
    et(i).e = [];
    et(i).n = 0;
    for j = 1:V
        if (A(i,j) ~= 0)
            et(i).e = [et(i).e j];
            et(i).n = et(i).n + 1;
        end
    end
    et(i).e=single(et(i).e);
end
end