clear all;
close all; 
clc;
NP=50; %������
Pc = 0.8; %ѡ�����
Pm = 0.05; %�������
f=30000000;%�ź�Ƶ��
c=30000000; 
lamda = c/f;%����
d = lamda/2;%���
beta = 2*pi/lamda;%����
seta0 = 0;%��������
G = 100;  %����Ŵ�����
L = 128; %������
NL = 16; %��Ԫ��
NN = 1800; %���̶ֿ�
E0 = 1; %��ѹ����
%%%���ɳ�ʼ��Ⱥ
f = randint(NP,L);  %����һ��L*NP�����������̬�ֲ�
%��ʼ��Ⱥ�ķֲ�
maxE = 1;
for i = 1:NP
    for j = 1:L
        if f(i,j)==1
            plot(i,j)
           hold on
        end
    end
end
%%%%�Ŵ��㷨ѭ��
for k = 1:G
    k
    %���룬����ʵ���� 
    % �����ֵ�԰�ȣ�����Ӧ��
    for i = 1:NP
        F(i) = -MSLL3(d,beta,NN,NL,seta0,f(i,:),L,maxE); % ȡ�����ģ����淽��ȡ���ֵ���Ҵ���
    end% ȡ��NP��������԰��ƽ
    maxF = max(F)
    minF = min(F)
    rr = find(F == maxF);
    fBest = f(rr(1),:);
    Fin = (F-minF)/(maxF-minF); %��һ��
    %rr(1,1)
    %fBest
    %maxFit
    %%%�����ֶ��̵�ѡ�����
    sum_Fin = sum(Fin);
    fitvalue = Fin./sum_Fin;
    fitvalue = cumsum(fitvalue);%cumsum(A) ����һ����Aͬ��ͬ�еľ���
                                %�����е�m�е�n��Ԫ����A�е�1�е���m
                                %�е����е�n��Ԫ�ص��ۼӺ�
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
    %���ʽ���
    for i = 1:2:NP
        P=rand;
        if P<Pc
            q = randint(1,L);%����1*L�������0��1
            for j = 1:L
                if q(j) == 1
                    temp = nf(i+1,j);
                    nf(i+1,j) = nf(i,j);
                    nf(i,j) = temp;
                end
            end
        end
    end
    %%%���ʱ���
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
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')
save fBestn.mat fBest
save L.mat L
save E.mat E0


 
