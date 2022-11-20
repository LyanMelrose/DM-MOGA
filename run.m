
%%
clc;
clear ;
Notrue=0;       %������ʵ���֡�1����ʵ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/5/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge_add=1;   %�߱��Ƿ��1��ʼ���Ǵ�0��ʼ����0��ʼedge_add = 1
% change=0;       %�����ڽӾ�����ʽ���Ǳ߱���ʽ   --1Ϊ����-----0Ϊ�߱�-------

%% SECTION TITLE
% DESCRIPTIVE TEXT
k=1;           %k���Ʋ�ͬ�����磬0Ϊ��ʵ�����磬1Ϊ5000����Ĵ���������磬2Ϊ1000������������
switch k
    %6����ʵ����;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/5/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
        Date={'karate','FB50','SFI','blogs','jazz','CA-HepTh_9877_reIndex','ca-HepTh2_12008_reIndex','CQC','netscience','hepth','hepth1','4000','jazz','protein-protein','football','dolphin','125_point','polbooks'};
        for i=3:3
            path = sprintf('RealWorld/%s.txt',Date{i});
            real_path=sprintf('RealWorld/real_label_%s.txt',Date{i});
            name=Date{i};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/5/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MAIN_START(path,real_path,name,Notrue,k);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/5/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ʵ��ComData
    case 1
%         Notrue=1;       %������ʵ���֡�
        Date={'facebook','github','karate_datatype_test'};
        for i=1:1
            path = sprintf('ComData/%s.mat',Date{i});
            name=Date{i};
            MAIN_START(path,'',name,Notrue,k);
        end
    
    %% ʵ��5000_20_100_0.0_0.8_20_50
    case 2
        Date={'0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        %%ǰ������������0 0 0��
        for i=5:9
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('new_LFR_5000_20_100/5000_%s.txt',Date{i});
            real_path=sprintf('new_LFR_5000_20_100/real_5000_%s.txt',Date{i});
            name=sprintf('new_LFR_5000_20_100_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end
        
        %% ʵ��1000_10_50_20_100_0.1_0.8����һ�����һ�û�ġ���2020/5/23��
    case 3
        %%
        Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'};
        for i=4:9
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR_1000_20_50_20_100/1000_%s.txt',Date{i});
            real_path=sprintf('LFR_1000_20_50_20_100/real_1000_%s.txt',Date{i});
            name=sprintf('LFR_1000_20_50_20_100_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end

        %%
    case 4
        Date={'1000','2000','3000','4000','5000'};
        for i=5:5
            path=sprintf('LFR_1000-5000/%s.txt',Date{i});
            real_path=sprintf('LFR_1000-5000/real_%s.txt',Date{i});
            name=sprintf('LFR_1000-5000_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020/5/23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% ʵ��LFR_200
    case 5
        Date={'0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'};
        for i=1:11
            path=sprintf('LFR_200/%s.txt',Date{i});
            real_path=sprintf('LFR_200/real%s.txt',Date{i});
            name=sprintf('LFR_200_%s.txt',Date{i});
            MAIN_START(path,real_path,name,Notrue,k);
        end
        
   %% 115����ķ���
    case 6
        Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        for i=7:8
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR_115/115_%s.txt',Date{i});
            real_path=sprintf('LFR_115/real_115_%s.txt',Date{i});
            name=sprintf('LFR_115_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end

   %% 25000����ķ��������������
    case 7
        %     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        Date={'0.1'};
        for i=1:1
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR_EXTEND/25000/%s/network.dat',Date{i});
            real_path=sprintf('LFR_EXTEND/25000/%s/community.dat',Date{i});
            name=sprintf('LFR_EXTEND_25000_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end
        
        %% # Directed graph (each unordered pair of nodes is saved once): Email-Enron.txt
        %% # Enron email network (edge indicated that email was exchanged, undirected edges)
        %% # Nodes: 36692 Edges: 367662
    case 8
        Date={'p2p-Gnutella25_undirect'};
        for i=1:1
            if i ==1
                edge_add =1;%%��0��ʼ
                Notrue=1;%%����ʵ����
            end
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('largeScale/%s.txt',Date{i});
            real_path=sprintf('');
            name=sprintf('largeScale_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end

    case 9
        Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        for i=5:5
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR_EXTEND/1000/1000_%s.txt',Date{i});
            real_path=sprintf('LFR_EXTEND/1000/1000_%s_labels.txt',Date{i});
            name=sprintf('LFR_1000_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end

        %% ʵ��1000_10_50_20_100_0.1_0.8
    case 1000
        Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        for i=6 :6
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR/1000/%s/network.dat',Date{i});
            real_path=sprintf('LFR/1000/%s/community.dat',Date{i});
            name=sprintf('1000_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end
        
    case 10
        Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
        for i=7:8
            clearvars -EXCEPT Date i change Noture edge_add;
            path=sprintf('LFR_EXTEND/5000/%s/network.dat',Date{i});
            real_path=sprintf('LFR_EXTEND/5000/%s/community.dat',Date{i});
            name=sprintf('LFR_EXTEND_5000_%s',Date{i});
            MAIN_START(path,name,real_path,change,Notrue,edge_add);
        end
end