% Input: 
%       (1) n: the number of nodes in G
%       (2) m: the size of WMAN, mKM x mKM
%       (3) up: the max distance between two APs
%       (4) low: the min distance between two APs
%       (5) mode: the computation demands of APs are subject to possion<0> or
%           pareto distribution<1>; default value is 0

% Output: the WMAN model denoted by model
%         model.size 表示节点个数
%         model.area 表示WMAN的区域大小
%         model.upBound 表示AP之间的距离上限
%         model.lowBound 表示AP之间的距离下限
%         model.distribution 表示AP的计算需求量的概率分布
%         model.locMat 表示AP的位置信息，即节点的坐标
%         model.adjMat 表示WMAN的邻接矩阵
%         model.G  表示将WMAN用图G表示
%         model.G.Nodes.demand  表示AP的计算需求量
%         model.hopMat 表示图的最短路径矩阵
%         model.capacity  表示边缘服务器的性能

% 输入检测->模型初始化->生成所需的70%的城市中心区域节点（点与点之间距离在[low,up]之间，故分布较为密集）
% ->生成所需的30%的郊区节点（点与点之间距离在[2*low,2*up]之间，故分布较为分散,注意应该让2*low>up）
% ->将上述生成节点保存到容器，中心区域距离在
function model = initModel(n,m,up,low,mode,capacity)
    
    % 输入合法性检测
    if up<low
        fprintf('up must be bigger than low');
        return;
    end
    % 模型初始化
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
    
    % 初始化n个节点的容器节点坐标，行为节点的个数，列为第i行的x坐标值和y坐标值
    Z = zeros(n,2);
    % 随机生成一个在m x m区域内的节点，添加至容器中，即第一个随机生成的节点
    % 为了统一，假设区域的中心为初始点
    T = [0.5,0.5]*m;
    Z(1,:) = T;
    
    
    % 生成节点的拓扑结构
    % 为了符合审稿人的要求， 生成的前百分之70的点之间的连接距离为后生成的百分之三十的点的一半
    % 生成前百分之70的节点，由于容器已经有了中心节点，节点索引从2开始
    for i = 2:(0.7*n)
        up_flag = true; 
        low_flag = true;
       
        % 产生符合要求的点
        while up_flag||low_flag
            %生成一个随机节点，其坐标为T
            T = rand(1,2)*m;
            %将新生成节点的距离上下界标识初始化，越上下界则标识为true，未越界标识为false，默认越上界
            up_flag = true; 
            low_flag = false;
            
            %将生成新节点的坐标T与Z中的已生成的i-1个节点坐标对比，需满足
            for j = 1:i-1
                 % 只要新生成的点距离现有的点中的某个点小于上界即可
                 if pdist([T;Z(j,:)],'euclidean') < up
                    up_flag = false;
                 end
                 % 新生成的点与所有的点的距离必须大于下界   
                 if pdist([T;Z(j,:)],'euclidean') < low
                    low_flag = true;
                 end
            end  
        end
        
        %跳出while循环则说明上下界表示皆为false，也就是未越界，这样我们就把T坐标加入到Z坐标的第i行，也就是第i个
        Z(i,:) = T;
        %disp(['运行到了',num2str(i),'了'])
    end
    
     %对剩余的百分之30的节点进行坐标生成，距离上下限为[2*low,2*up]
     for i = (0.7*n):n
        up_flag = true; 
        low_flag = true;
       
       
        % 产生符合要求的点
        while up_flag||low_flag
            T = rand(1,2)*m;
            up_flag = true; 
            low_flag = false;
            
            for j = 1:i-1
                 % 只要新生成的点距离现有的点中的某个点小于上界即可
                 if  pdist([T;Z(j,:)],'euclidean') < 2*up 
                    up_flag = false;
                 end
                 % 新生成的点与所有的点的距离必须大于下界   
                 if pdist([T;Z(j,:)],'euclidean') < 2*low 
                    low_flag = true;
                 end
            end
            
             
        end
%         disp(['运行到了',num2str(i),'了'])
        Z(i,:) = T;
        
    end
    % 生成节点拓扑结构结束
    
    % 将节点的位置信息保存至模型的参数locMat中
    model.locMat = Z;
    
    % 节点之间的距离信息
    model.disMat = squareform(pdist(Z));
    % 邻接矩阵信息
    model.adjMat = (model.disMat < up) - eye(n);
    model.adjMat = zeros(n);
    model.adjMat(1:(0.7*n),1:(0.7*n)) = (model.disMat(1:(0.7*n),1:(0.7*n)) < up);
    model.adjMat(1:(0.7*n),(0.7*n)+1:n) = (model.disMat(1:(0.7*n),(0.7*n)+1:n) < 2*up);
    model.adjMat((0.7*n)+1:n,(0.7*n)+1:n) = (model.disMat((0.7*n)+1:n,(0.7*n)+1:n) < 2*up);
    model.adjMat((0.7*n)+1:n,1:(0.7*n)) = (model.disMat((0.7*n)+1:n,1:(0.7*n)) < 2*up);
    
    % 生成每个AP的计算需求量
    % 此处暂时不考虑泊松分布抑或帕累托分布，统一使用XU等人使用的计算量生成方法
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
     % 最短路径矩阵信息
    model.hopMat = distances(model.G);
end


