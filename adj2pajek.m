% Convert a symmetrical adjacency matrix to a Pajek network file
%
% Input
%   - adj : adjacency matrix
%   - community : 聚类后的结果
%   - community_num : 想要画的社团的的编号,它是一个数组
%   - name: file name (without extension)
%   - dir : destination directory
%
% Author: Erwan Le Martelot
% Date: 31/05/11

function [] = adj2pajek(adj , community , community_num , name, dir)
% adj=[ 0     1     1     1     0     0     1
%      1     0     0     1     1     0     0
%      1     0     0     1     1     0     1
%      1     1     1     0     1     1     0
%      0     1     1     1     0     1     0
%      0     0     0     1     1     0     0
%      1     0     1     0     0     0     0];
% 
%  community={[1 2 3 4 5 6 7]};
%  community_num=1;
%  name='clique_S3';
%  dir='C:\Documents and Settings\Administrator\桌面\网络画图';
 

    dst_file = [name,'.net'];
    nargin=5;
    if nargin == 5
        dst_file = [dir, '\', dst_file];
    end
    
    fid = fopen(dst_file, 'w');                               %文件路径打开不能打开已经链接的
    
    %如果community_num为all，则默认为画出所有的聚类结果，否则是只画community中的一小部分
    if strcmp(community_num,'all')
        
        fprintf(fid, '*Vertices %d\n', size(adj,1));

    %     for i=1:size(adj,1)
    %         fprintf(fid, ' %d "%d"\n', i, vertices_num{i});
    %     end
    
%         a=[];
%         for i=1:length(community_num)
%             a = union(a,community{community_num(i)});
%         end;
        
        for i=1:size(adj,1)
            fprintf(fid, ' %d %d\n', i, i);
        end

        fprintf(fid, '*Arcs\n');
        fprintf(fid, '*Edges\n');

        for i=1:size(adj,1)
            for j=i+1:size(adj,1)
                if adj(i,j) > 0
                    fprintf(fid, ' %d %d %f\n', i, j, adj(i,j));
                end
            end
        end
    else
        a=[];
        for i=1:length(community_num)
            a = union(a,community{community_num(i)});
        end
        
%         %去掉度数为1的点
%         i=1;
%         while i<=length(a)
%             if degree(a(i))==1
%                 a(i)=[];
%                 i=i-1;
%             end;
%             i=i+1;
%         end;

        fprintf(fid, '*Vertices %d\n', length(a));

    %     for i=1:size(adj,1)
    %         fprintf(fid, ' %d "%s"\n', i, vertices_num{i});
    %     end

        for i=1:length(a)
            fprintf(fid, '%d "%d"\n', i, i);
        end

        fprintf(fid, '*Arcs\n');
        fprintf(fid, '*Edges\n');

        for i=1:length(a)
            for j=i+1:length(a)
                m=a(i)
                n=a(j)
                if adj(m,n)==1
%                 if ismember(a(i),adj(a(j)+1))               %%判断是否是子集
                    fprintf(fid, '%d %d %f\n', i, j, 1);
                end
            end
        end
    end

    fclose(fid);