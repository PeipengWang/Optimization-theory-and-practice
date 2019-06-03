clear all;
close all; 
clc;
load('fBestn.mat');
load('L.mat')
load('E.mat')
fbest = 0;
NP=50; %������
Pc = 0.8; %ѡ�����
Pm = 0.01; %�������
f=30000000;%�ź�Ƶ��
c=30000000; 
lamda = c/f;%����
d = lamda/2;%���
beta = 2*pi/lamda;%����
seta0 = 0;%��������
G = 200;  %����Ŵ�����
NL = 16; %��Ԫ��
NN = 1800; %���̶ֿ�
e = 1;
maxE = 8;
seta = linspace(-pi,pi,NN);
%for m = 1:NN
    %alfa = sum(2*pi*fBest*(2.^(0:L-1)')/(2^L-1)); %����
    %fai = alfa+beta*d.*cos(seta(m));
    %F1(m) = abs(sin(NL*(fai./2))./sin(fai./2));
%end
alfa = (2*pi*(fBest(1:L/2)*(2.^(0:L/2-1)'))/(2^(L/2)-1)); %����
x = E0*fBest(L/2+1:L)*(2.^(0:L/2-1)')/(2^(L/2)-1); %����

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

% ��һ������ͼ
FdB = FdB-3;
figure
plot(seta*180/pi,FdB)
xlabel('\theta/(��)')
ylabel('��������/dB')
axis([0 180 -50 0]);%��Ϊx���y������һ�����޷�Χ
grid on
%%%����󸱰��ƽ

for i = 1:2
    maxFdB = max(FdB);
    rr = find(FdB == maxFdB); %find()�����Ĺ������ҵ��������߾����в�Ϊ0��Ԫ�أ�
    maxFdB = max(FdB);
    rr = find(FdB == maxFdB); %find()�����Ĺ������ҵ��������߾����в�Ϊ0��Ԫ�أ�
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
MSLL = max(FdB)                      %�������Ҫ�ҵ���������һ��������Ԫ������
figure 
polar(seta,FdB)
xlabel('\theta/(��)')
ylabel('��������/dB')


