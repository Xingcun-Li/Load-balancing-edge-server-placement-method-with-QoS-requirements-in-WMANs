% Input: 
%       (1) n: the number of nodes in G
%       (2) m: the size of WMAN, mKM x mKM
%       (3) up: the max distance between two APs
%       (4) low: the min distance between two APs
%       (5) mode: the computation demands of APs are subject to possion<0> or
%           pareto distribution<1>; default value is 0

% Output: the WMAN model denoted by model
%         model.size ��ʾ�ڵ����
%         model.area ��ʾWMAN�������С
%         model.upBound ��ʾAP֮��ľ�������
%         model.lowBound ��ʾAP֮��ľ�������
%         model.distribution ��ʾAP�ļ����������ĸ��ʷֲ�
%         model.locMat ��ʾAP��λ����Ϣ�����ڵ������
%         model.adjMat ��ʾWMAN���ڽӾ���
%         model.G  ��ʾ��WMAN��ͼG��ʾ
%         model.G.Nodes.demand  ��ʾAP�ļ���������
%         model.hopMat ��ʾͼ�����·������
%         model.capacity  ��ʾ��Ե������������

% ������->ģ�ͳ�ʼ��->���������70%�ĳ�����������ڵ㣨�����֮�������[low,up]֮�䣬�ʷֲ���Ϊ�ܼ���
% ->���������30%�Ľ����ڵ㣨�����֮�������[2*low,2*up]֮�䣬�ʷֲ���Ϊ��ɢ,ע��Ӧ����2*low>up��
% ->���������ɽڵ㱣�浽�������������������
function model = initModel(n,m,up,low,mode,capacity)
    
    % ����Ϸ��Լ��
    if up<low
        fprintf('up must be bigger than low');
        return;
    end
    % ģ�ͳ�ʼ��
    model.size = n;
    model.rawNodes = 1:n;
    model.area = m; 
    model.upBound = up;
    model.lowBound = low;
    model.distribution = 'possion';
    if mode==1
    	model.distribution = 'pareto';
    end
    model.capacity = capacity;
    
    % ��ʼ��n���ڵ�������ڵ����꣬��Ϊ�ڵ�ĸ�������Ϊ��i�е�x����ֵ��y����ֵ
    Z = zeros(n,2);
    % �������һ����m x m�����ڵĽڵ㣬����������У�����һ��������ɵĽڵ�
    % Ϊ��ͳһ���������������Ϊ��ʼ��
    T = [0.5,0.5]*m;
    Z(1,:) = T;
    
    
    % ���ɽڵ�����˽ṹ
    % Ϊ�˷�������˵�Ҫ�� ���ɵ�ǰ�ٷ�֮70�ĵ�֮������Ӿ���Ϊ�����ɵİٷ�֮��ʮ�ĵ��һ��
    % ����ǰ�ٷ�֮70�Ľڵ㣬���������Ѿ��������Ľڵ㣬�ڵ�������2��ʼ
    for i = 2:(0.7*n)
        up_flag = true; 
        low_flag = true;
       
        % ��������Ҫ��ĵ�
        while up_flag||low_flag
            %����һ������ڵ㣬������ΪT
            T = rand(1,2)*m;
            %�������ɽڵ�ľ������½��ʶ��ʼ����Խ���½����ʶΪtrue��δԽ���ʶΪfalse��Ĭ��Խ�Ͻ�
            up_flag = true; 
            low_flag = false;
            
            %�������½ڵ������T��Z�е������ɵ�i-1���ڵ�����Աȣ�������
            for j = 1:i-1
                 % ֻҪ�����ɵĵ�������еĵ��е�ĳ����С���Ͻ缴��
                 if pdist([T;Z(j,:)],'euclidean') < up
                    up_flag = false;
                 end
                 % �����ɵĵ������еĵ�ľ����������½�   
                 if pdist([T;Z(j,:)],'euclidean') < low
                    low_flag = true;
                 end
            end  
        end
        
        %����whileѭ����˵�����½��ʾ��Ϊfalse��Ҳ����δԽ�磬�������ǾͰ�T������뵽Z����ĵ�i�У�Ҳ���ǵ�i��
        Z(i,:) = T;
        %disp(['���е���',num2str(i),'��'])
    end
    
     %��ʣ��İٷ�֮30�Ľڵ�����������ɣ�����������Ϊ[2*low,2*up]
     for i = (0.7*n):n
        up_flag = true; 
        low_flag = true;
       
       
        % ��������Ҫ��ĵ�
        while up_flag||low_flag
            T = rand(1,2)*m;
            up_flag = true; 
            low_flag = false;
            
            for j = 1:i-1
                 % ֻҪ�����ɵĵ�������еĵ��е�ĳ����С���Ͻ缴��
                 if  pdist([T;Z(j,:)],'euclidean') < 2*up 
                    up_flag = false;
                 end
                 % �����ɵĵ������еĵ�ľ����������½�   
                 if pdist([T;Z(j,:)],'euclidean') < 2*low 
                    low_flag = true;
                 end
            end
            
             
        end
%         disp(['���е���',num2str(i),'��'])
        Z(i,:) = T;
        
    end
    % ���ɽڵ����˽ṹ����
    
    % ���ڵ��λ����Ϣ������ģ�͵Ĳ���locMat��
    model.locMat = Z;
    
    % �ڵ�֮��ľ�����Ϣ
    model.disMat = squareform(pdist(Z));
    % �ڽӾ�����Ϣ
    model.adjMat = (model.disMat < up) - eye(n);
    model.adjMat = zeros(n);
    model.adjMat(1:(0.7*n),1:(0.7*n)) = (model.disMat(1:(0.7*n),1:(0.7*n)) < up);
    model.adjMat(1:(0.7*n),(0.7*n)+1:n) = (model.disMat(1:(0.7*n),(0.7*n)+1:n) < 2*up);
    model.adjMat((0.7*n)+1:n,(0.7*n)+1:n) = (model.disMat((0.7*n)+1:n,(0.7*n)+1:n) < 2*up);
    model.adjMat((0.7*n)+1:n,1:(0.7*n)) = (model.disMat((0.7*n)+1:n,1:(0.7*n)) < 2*up);
    
    % ����ÿ��AP�ļ���������
    % �˴���ʱ�����ǲ��ɷֲ��ֻ������зֲ���ͳһʹ��XU����ʹ�õļ��������ɷ���
    if mode == 0
%           model.comMat =  poissrnd(5,1,n);
          num = randi([50,500],n,1);
          comUnit = randi([50,200],n,1);
          model.comMat = num.*comUnit;
    else 
%         model.comMat =  gprnd(1,2,3,1,n);
        num = randi([50,500],n,1);
        comUnit = randi([50,200],n,1);
        model.comMat = num.*comUnit;
    end
    model.G = graph(model.adjMat);
    model.G.Nodes.demand = model.comMat; 
     % ���·��������Ϣ
    model.hopMat = distances(model.G);
end


