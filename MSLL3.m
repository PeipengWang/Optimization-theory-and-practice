%%%��������԰�%%%
function  MSLL3 = MSLL3(d,beta,NN,NL,seta0,f,L,maxE)
seta = linspace(-pi,pi,NN);
alfa = (2*pi*(f(1:L/2)*(2.^(0:L/2-1)'))/(2^(L/2)-1)); %����
x = maxE/8*f(L/2+1:L)*(2.^(0:L/2-1)')/(2^(L/2)-1); %����

I1 = linspace(maxE-(NL/2)*x,maxE,NL/2);
I2 = linspace(maxE,maxE-(NL/2*x),NL/2);
I = [I1,I2]-x;
for m = 1:NN
    fai = beta*d*([0:NL-1]+1-(NL+1)/2)*cos(seta(m)-seta0)-[0:NL-1]*alfa;
    Fit(m) = abs(sum(I*exp(sqrt(-1)*(fai)')));
end
%plot(alfa,fit)
%hold on
FdB = 20*log10(Fit/max(Fit));

%�������Ҫ�ҵ���������һ��������Ԫ������
%%%ȡ�԰��ƽ�������ֵ֮��Ĵδ�ֵ
FdB = FdB-3;
for i = 1:2
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
MSLL3 = max(FdB);
end