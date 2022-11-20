% Convert a community list to a Pajek partition file
%
% Input
%   - com : community list per vertex
%   - community_num : 想要画出来的社团编号，它是一个数组
%   - name: file name (without extension)
%   - dir : destination directory
%
% Author: Erwan Le Martelot
% Date: 01/06/11


function [] = com2pajek(community , community_num , name, dir)

    dst_file = [name,'.clu'];
    if nargin == 4
        dst_file = [dir, '\', dst_file];
    end

    fid = fopen(dst_file, 'w');
    
    if strcmp(community_num,'all')
        a=[];
        for i=1:length(community)
            a=union(a,community{i});
        end;
        vertex_num=length(a);

        fprintf(fid, '*Vertices %d\n', vertex_num);

        x=[];
        b=length(community)+1;
         for node=1:vertex_num
             a1=0;
             time=0;
            for k=1:length(community)
                com2=community{k};
                if ismember(a(node),com2)
                    a1=k;
                    time=time+1;
                end;
            end;

            if time==1
                x=[x;a1];
            else
                x=[x;b];
            end;
        end;

        for i=1:vertex_num
            fprintf(fid,'%d\n',x(i));
        end;
    else
        a=[];
        for i=1:length(community_num)
            a=union(a,community{community_num(i)});
        end;
        
%         %去掉度数为1的点
%         i=1;
%         while i<=length(a)
%             if degree(a(i))==1
%                 a(i)=[];
%                 i=i-1;
%             end;
%             i=i+1;
%         end;
        vertex_num=length(a);

        fprintf(fid, '*Vertices %d\n', vertex_num);

        c=[];
        b=length(community_num)+1;
         for node=1:vertex_num
             a1=0;
             time=0;
            for k=1:length(community_num)
                com2=community{community_num(k)};
                if ismember(a(node),com2)
                    a1=k;
                    time=time+1;
                end;
            end;

            if time==1
                c=[c;a1];
            else
                c=[c;b];
            end;
        end;

        for i=1:vertex_num
            fprintf(fid,'%d\n',c(i));
        end;
    end;
    
    fclose(fid);

end
