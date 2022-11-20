clear 
clc
tic
global numObjectives popsize numVar niche max_gen mutate_posibility  edgeslist V  %��Ŀ��ȫ�ֲ���
global strNetwork AdjacentMatrix Datalabel  idealp weights neighbors  ll M  adj_mat  pc pg pm             %����������
global Thirty_Run_maxQ Thirty_Run_maxNMI      %ʵ�����
%% ��Ŀ����������ֵ

max_gen         = 100;       % maximal number of generations
niche           = 40;        % neighborhood size
popsize			= 100;		 % population size
runtimes        = 10;		 % run times
M = 2;
pg=1;
pm=0.9;
pc=1;

%% ������� numVar;�ڽӾ�����ʵ���֣�ÿ�������Ϣnode

strNetwork='karate';
AdjacentMatrix = load('RealWorld\karate.txt');
Datalabel=load('RealWorld\real_label_karate.txt');
adj_mat=AdjacentMatrix;
[V,V] = size(adj_mat); 
edgeslist = edges_list(adj_mat,V);

numVar=size(AdjacentMatrix,1);
 ll=0;
for i=1:V
    for j=i:V
        ll=adj_mat(i,j)+ll;
    end
end

%% ʵ�鲿��
root = sprintf('results1/%s/ParetoFront',strNetwork);  %%�����ļ�����ʵ������
if ~isdir(root) %�ж�·���Ƿ����
    mkdir(root); 
end
root = sprintf('results1/%s/metrics',strNetwork);
if ~isdir(root) %�ж�·���Ƿ����
    mkdir(root);
end

Thirty_Run_maxQ=[];
Thirty_Run_maxNMI=[];

for mg=1:runtimes
    fprintf('��%s������\n',num2str(mg));
    TMOEAD(mg);
    fprintf('\n');
end 

%[~,index1]=max(Thirty_Run_maxQ(:,1));%�ҵ�����Q����
%[~,index2]=max(Thirty_Run_maxNMI(:,2));%�ҵ�����NMI����

fprintf('Thirtyrun_results:\n');
fprintf('maxQ =    %g %g\n',Thirty_Run_maxQ(index1,1),Thirty_Run_maxQ(index1,2));    %�������Qʱ��Q��NMI
fprintf('maxNMI =  %g %g\n',Thirty_Run_maxNMI(index2,1),Thirty_Run_maxNMI(index2,2));%�������NMIʱ��NMI��Q

path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxQ.txt',strNetwork,strNetwork);%����ʵ������ȡƽ��ֵ
savedata1(path,[Thirty_Run_maxQ;0 0;mean(Thirty_Run_maxQ(:,1)) mean(Thirty_Run_maxQ(:,2))]); 
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxNMI.txt',strNetwork,strNetwork);
savedata1(path,[Thirty_Run_maxNMI;0 0;mean(Thirty_Run_maxNMI(:,1)) mean(Thirty_Run_maxNMI(:,2))]); 
