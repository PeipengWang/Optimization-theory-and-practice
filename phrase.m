clear all;
close all; 
clc;
NP=50; %个体数
Pc = 0.8; %选择概率
Pm = 0.05; %变异概率
f=30000000;%信号频率
c=30000000; 
lamda = c/f;%波长
d = lamda/2;%间距
beta = 2*pi/lamda;%波数
seta0 = 0;%波束方向
G = 100;  %最大遗传代数
L = 128; %编码数
NL = 16; %阵元数
NN = 1800; %划分刻度
E0 = 1; %电压电流
%%%生成初始种群
f = randint(NP,L);  %产生一个L*NP的随机矩阵，正态分布
%初始种群的分布
maxE = 1;
for i = 1:NP
    for j = 1:L
        if f(i,j)==1
            plot(i,j)
           hold on
        end
    end
end
%%%%遗传算法循环
for k = 1:G
    k
    %解码，个体实际数 
    % 计算峰值旁瓣比，即适应度
    for i = 1:NP
        F(i) = -MSLL3(d,beta,NN,NL,seta0,f(i,:),L,maxE); % 取成正的，下面方便取最大值并且处理
    end% 取得NP个个体的旁瓣电平
    maxF = max(F)
    minF = min(F)
    rr = find(F == maxF);
    fBest = f(rr(1),:);
    Fin = (F-minF)/(maxF-minF); %归一化
    %rr(1,1)
    %fBest
    %maxFit
    %%%基于轮赌盘的选择操作
    sum_Fin = sum(Fin);
    fitvalue = Fin./sum_Fin;
    fitvalue = cumsum(fitvalue);%cumsum(A) 返回一个和A同行同列的矩阵，
                                %矩阵中第m行第n列元素是A中第1行到第m
                                %行的所有第n列元素的累加和
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi <= NP
        if(ms(newi)<fitvalue(fiti))
            nf(newi,:) = f(fiti,:);
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end
    %概率交叉
    for i = 1:2:NP
        P=rand;
        if P<Pc
            q = randint(1,L);%生成1*L随机矩阵，0或1
            for j = 1:L
                if q(j) == 1
                    temp = nf(i+1,j);
                    nf(i+1,j) = nf(i,j);
                    nf(i,j) = temp;
                end
            end
        end
    end
    %%%概率变异
    for m = 1:NP
        for n = 1:L
            r = rand(1,1);
            if r < Pm
                nf(m,n) = ~nf(m,n);
            end
        end
    end 
    f = nf;
    f(1,:) = fBest;
    trace(k) = maxF;

end
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')
save fBestn.mat fBest
save L.mat L
save E.mat E0


 
