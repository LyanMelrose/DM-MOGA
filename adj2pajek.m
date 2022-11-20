% Convert a symmetrical adjacency matrix to a Pajek network file
%
% Input
%   - adj : adjacency matrix
%   - community : �����Ľ��
%   - community_num : ��Ҫ�������ŵĵı��,����һ������
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
%  dir='C:\Documents and Settings\Administrator\����\���续ͼ';
 

    dst_file = [name,'.net'];
    nargin=5;
    if nargin == 5
        dst_file = [dir, '\', dst_file];
    end
    
    fid = fopen(dst_file, 'w');                               %�ļ�·���򿪲��ܴ��Ѿ����ӵ�
    
    %���community_numΪall����Ĭ��Ϊ�������еľ�������������ֻ��community�е�һС����
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
        
%         %ȥ������Ϊ1�ĵ�
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
%                 if ismember(a(i),adj(a(j)+1))               %%�ж��Ƿ����Ӽ�
                    fprintf(fid, '%d %d %f\n', i, j, 1);
                end
            end
        end
    end

    fclose(fid);