clear all;
close all; 
clc;
load('fBestn.mat');
load('L.mat')
load('E.mat')
fbest = 0;
NP=50; %个体数
Pc = 0.8; %选择概率
Pm = 0.01; %变异概率
f=30000000;%信号频率
c=30000000; 
lamda = c/f;%波长
d = lamda/2;%间距
beta = 2*pi/lamda;%波数
seta0 = 0;%波束方向
G = 200;  %最大遗传代数
NL = 16; %阵元数
NN = 1800; %划分刻度
e = 1;
maxE = 8;
seta = linspace(-pi,pi,NN);
%for m = 1:NN
    %alfa = sum(2*pi*fBest*(2.^(0:L-1)')/(2^L-1)); %解码
    %fai = alfa+beta*d.*cos(seta(m));
    %F1(m) = abs(sin(NL*(fai./2))./sin(fai./2));
%end
alfa = (2*pi*(fBest(1:L/2)*(2.^(0:L/2-1)'))/(2^(L/2)-1)); %解码
x = E0*fBest(L/2+1:L)*(2.^(0:L/2-1)')/(2^(L/2)-1); %解码

I1 = linspace(maxE-(NL/2)*x,maxE,NL/2);
I2 = linspace(maxE,maxE-(NL/2*x),NL/2);
I = [I1,I2]-x;
for m = 1:NN
    fai = beta*d*([0:NL-1]+1-(NL+1)/2)*cos(seta(m)-seta0)-[0:NL-1]*alfa;
    F1(m) = abs(sum(I*exp(sqrt(-1)*(fai)')));
end
F1 = abs(F1);
D1 = 0; 
for m = 1:NL-1;
    D1 = (NL-m)/(m*beta*d)*sin(m*beta*d)*cos(m*alfa); 
end
D = 1/(1/NL+(2/NL^2)*D1)


%plot(alfa,fit)
%hold on

plot(seta,F1);
FdB = 20*log10(F1/max(F1));

% 归一化方向图
FdB = FdB-3;
figure
plot(seta*180/pi,FdB)
xlabel('\theta/(度)')
ylabel('阵列增益/dB')
axis([0 180 -50 0]);%可为x轴和y轴设置一个极限范围
grid on
%%%求最大副瓣电平

for i = 1:2
    maxFdB = max(FdB);
    rr = find(FdB == maxFdB); %find()函数的功能是找到向量或者矩阵中不为0的元素，
    maxFdB = max(FdB);
    rr = find(FdB == maxFdB); %find()函数的功能是找到向量或者矩阵中不为0的元素，
    mm = rr(1);
    tu_up = 0;
    while(FdB(mm+tu_up) >= FdB(mm+tu_up+1))
        tu_up = tu_up+1;
        if mm+tu_up+1>1800
            mm = -tu_up+1;
        end
    end
    mm = rr(1);
    tu_down = 0;
    while (FdB(mm-tu_down) >= FdB(mm-tu_down-1))
        tu_down = tu_down+1;
        if mm-tu_down-1<1
            mm  = 1800+tu_down;
        end
    end
    FdB(mm-tu_down:mm+tu_up) = -50;
end

figure(3)
plot(seta*180/pi,FdB);
axis([-180,180,-30,10])
MSLL = max(FdB)                      %那如果需要找到其中满足一定条件的元素索引
figure 
polar(seta,FdB)
xlabel('\theta/(度)')
ylabel('阵列增益/dB')


