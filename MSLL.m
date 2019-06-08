%%%��������԰�%%%
function  MSLL3 = MSLL3(d,beta,NN,NL,seta0,f,L,maxE)
seta = linspace(-pi,pi,NN);
%alfa = (2*pi*(f(1:L/2)*(2.^(0:L/2-1)'))/(2^(L/2)-1)); %����
alfa = 0;
%x = maxE/8*f(L/2+1:L)*(2.^(0:L/2-1)')/(2^(L/2)-1); %����
%I = [0.2035,0.2792,0.4089,0.5625,0.7147,0.8474,0.9466,1,1,0.9466,0.8474,...
%   0.7147,0.5625,0.4089,0.2792,0.2035];
for  i = 0:NL/2-1
    I1(i+1) = (maxE*(f(i*L/16+1:(i+1)*L/16)*(2.^(0:(L/16-1))'))/(2^(L/16)-1));
    %alfa1(i+1) = (2*pi*(f((i+8)*L/16+1:(i+9)*L/16)*(2.^(0:(L/16-1))'))/(2^(L/16)-1));
end
I2 = fliplr(I1);
I = [I1,I2];
for m = 1:NN
    fai = beta*d*([0:NL-1]+1-(NL+1)/2)*cos(seta(m)-seta0)-[0:NL-1]*alfa;
    Fit(m) = abs(sum(I*exp(sqrt(-1)*(fai)')));
end
%plot(alfa,fit)
%hold on
FdB = 20*log10(Fit/max(Fit));
%�������Ҫ�ҵ���������һ��������Ԫ������
%%%ȡ�԰��ƽ�������ֵ֮��Ĵδ�ֵ
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