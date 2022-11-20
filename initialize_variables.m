function f= initialize_variables(N,el,M,degree1,matrix,Node3)
% function [f]= initialize_variables(N,el,M,ll,degree1,matrix,Node3)
% �����matrix����ͨ�ڽӾ���������ľ���δ��
global idealp;
V=size(matrix,1);

K = M + V;

% NΪ�������������ÿ�����壬�����ھ���f�С�
% ÿ�д���һ�����ӣ�ÿ�д���һ���ھӾֲ�ģ���ĳ���ڵ㣬Ԫ���д洢��Ӧ�ڵ���ھӵı��
for i = 1 : N
    if rand<0.2
        f(i,1:V)=IGLP(V,el,Node3);
    else
        % elΪedgeslist
        for j = 1 : V
            % �ӵ�j���ֲ�ģ�������ѡ��һ���ڵ���ھӱ��
            f(i,j) = el(j).e(ceil(rand*el(j).n));
        end
    end
    %ÿ�е�������д洢Ŀ�꺯��ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2021/8/31%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(i,V + 1: K) = evaluate_objective(f(i,1:V),degree1,matrix,Node3);
%     f(i,V + 1: K) = evaluate_objective(f(i,1:V),ll,degree1,matrix,Node3);

end

% ��ʼ��ΪȺ�������ֱ���С��Ŀ�꺯��ֵ
idealp = min(f(:,V+1:K));
 
end
