% function [Sxy,alphao,fo] = crossfam(x,y,fs,df,dalpha)
fs=1280;
fc=160;
N=8192;
PW=N/fs;
signal_type1=2;
switch signal_type1
    case 1; input_signal = sin(2*pi*(fc)*((0:(N-1))/fs));
    case 2; input_signal = bpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 3; input_signal = qpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 4; input_signal = lfm_generator(fs);
    otherwise;
end
x = awgn(input_signal,10);
signal_type2=3;
switch signal_type2
    case 1; input_signal = sin(2*pi*(fc)*((0:(N-1))/fs));
    case 2; input_signal = bpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 3; input_signal = qpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 4; input_signal = lfm_generator(fs);
    otherwise;
end
y=awgn(input_signal,10);
M=128;
dalpha=fs/N;
df=dalpha*M;
Np = pow2(nextpow2(fs/df));     %输入信道数
L = Np/4;                       %同一行连续列的相邻点的偏移
P = pow2(nextpow2(fs/dalpha/L));%信道矩阵的列数
N = P*L;                        %输入数据的点数
if length(x) < N
    x(N) = 0;
elseif length(x) > N
    x = x(1:N);
end
if length(y) < N
    y(N) = 0;
elseif length(y) > N
    y = y(1:N);
end
NN = (P-1)*L+Np;
xx = x;
yy = y;
xx(NN) = 0;
yy(NN) = 0;
xx = xx(:);
yy = yy(:);
X = zeros(Np,P);
Y = zeros(Np,P);
for k = 0:P-1
    X(:,k+1) = xx(k*L+1:k*L+Np);
    Y(:,k+1) = yy(k*L+1:k*L+Np);
    %输入数据xx的1到Np存入X第一列，L+1到L+Np存入第二列，2L+1到2L+Np存入第三列，以此类推。即L为偏移，Np为每段长度。
end
a = hamming(Np);
XW = diag(a)*X;   %每段加窗
YW = diag(a)*Y;
%----------------First FFT--------------------%
XF1 = fft(XW);
YF1 = fft(YW);
% clear XW;   
XF1 = fftshift(XF1); 
YF1 = fftshift(YF1);
XF1 = [XF1(:,P/2+1:P) XF1(:,1:P/2)];  %这两行的功能是把傅里叶变换的结果上半块和下半块交换,左半块和右半块互换
YF1 = [YF1(:,P/2+1:P) YF1(:,1:P/2)];

%-------------Downconversion------------------%
E = zeros(Np,P);
for k = -Np/2+1:Np/2
    for m = 0:P-1
        E(k+Np/2,m+1) = exp(-i*2*pi*k*m*L/Np);
    end
end
XD = XF1.*E; 
YD = YF1.*E;
% clear XF1; 
XD = conj(XD'); %XD转置的复共轭
YD = conj(YD');

% clear ('XF1', 'E', 'XW', 'X', 'x'); 
%-----------Multiplication--------------------%
XYM = zeros(P,Np^2);
for k = 1:Np
    for l = 1:Np
        XYM(:,(k-1)*Np+l) = (XD(:,k).*conj(YD(:,l)));
    end
end
% clear XD;

%------------Second FFT-----------------------%
XYF2 = fft(XYM);  
% clear XM;
XYF2 = fftshift(XYF2);
XYF2 = [XYF2(:,Np^2/2+1:Np^2) XYF2(:,1:Np^2/2)];
% length(XF2);
XYF2 = XYF2(P/4:3*P/4,:);
M = abs(XYF2);
% clear XF2;  %Absolute value and complex magnitude

%%%%%%%%%%%%%%%%频率分辨率和循环频率分辨率
alphao = (-1:1/N:1)*fs;     %about the variable N!!!
fo = (-0.5:1/Np:0.5)*fs;
Sxy = zeros(Np+1,2*N+1);     %about the variable N!!!
for k1 = 1:P/2+1
    for k2 = 1:Np^2
        if rem(k2,Np) == 0
            l = Np/2;    
        else
            l = rem(k2,Np)-Np/2; 
        end
        k = ceil(k2/Np)-Np/2; 
        p = k1-P/4-1;
        alpha = (k-l)/Np+(p-1)/L/P;
        f = (k+l)/2/Np;
        if alpha < -1 | alpha > 1
            k2 = k2+1;
        elseif f < -0.5 | f > 0.5
            k2 = k2+1;
%         elseif rem(k+l,2)==0 && rem(1+N*(alpha+1),1)==0
        else
            kk = 1+Np*(f+0.5);
            ll = 1+N*(alpha + 1);
%           Sx(round(kk), round(ll)) = M(k1,k2);
            Sxy(round(kk), round(ll)) = M(k1,k2); 
        end
    end
end
Sxy = Sxy./max(max(Sxy)); % Normalizes the magnitudes of the values in output matrix (maximum = 1)
figure(1)
mesh(alphao,fo,Sxy)
title('SCD estimate using FAM')
xlabel('alpha');
ylabel('f');
zlabel('Sx');
figure(2);
plot(alphao,Sxy)
xlabel('Cycle frequency(alpha)')
ylabel('Sx')
figure(3);
plot(fo,Sxy)
xlabel('frequency(f)')
ylabel('Sx')
% figure(3)
% plot(fo,Sxy(:,1+N*((alpha/fs+1))));
% xlabel('f(Hz)');
% ylabel('Sx(alpha)');