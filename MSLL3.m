%%%计算最大旁瓣%%%
function  MSLL3 = MSLL3(d,beta,NN,NL,seta0,f,L,maxE)
seta = linspace(-pi,pi,NN);
alfa = (2*pi*(f(1:L/2)*(2.^(0:L/2-1)'))/(2^(L/2)-1)); %解码
x = maxE/8*f(L/2+1:L)*(2.^(0:L/2-1)')/(2^(L/2)-1); %解码

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

%那如果需要找到其中满足一定条件的元素索引
%%%取旁瓣电平，即最大值之外的次大值
FdB = FdB-3;
for i = 1:2
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
MSLL3 = max(FdB);
end