% function [SxY,alphao,fo] = crossssca(x,y,fs,df,dalpha)
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
Np = pow2(nextpow2(fs/df));     %�����ŵ���
L = Np/4;                       %ͬһ�������е����ڵ��ƫ��
P = pow2(nextpow2(fs/dalpha/L));%�ŵ����������
N = P*L;                        %�������ݵĵ���
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
xx(NN) = 0;
xx = xx(:);
X = zeros(Np,P);
for k = 0:P-1
    X(:,k+1) = xx(k*L+1:k*L+Np);
    %��������xx��1��Np����X��һ�У�L+1��L+Np����ڶ��У�2L+1��2L+Np��������У��Դ����ơ���LΪƫ�ƣ�NpΪÿ�γ��ȡ�
end
a = hamming(Np);
XW = diag(a)*X;   %ÿ�μӴ�
%----------------First FFT--------------------%
XF1 = fft(XW);
% clear XW;   
XF1 = fftshift(XF1); 
XF1 = [XF1(:,P/2+1:P) XF1(:,1:P/2)];  %�����еĹ����ǰѸ���Ҷ�任�Ľ���ϰ����°�齻��,������Ұ�黥��

%-------------Downconversion------------------%
E = zeros(Np,P);
for k = -Np/2+1:Np/2
    for m = 0:P-1
        E(k+Np/2,m+1) = exp(-i*2*pi*k*m*L/Np);
    end
end
XD = XF1.*E; 
% clear XF1; 
XD = conj(XD'); %XDת�õĸ�����;
% clear ('XF1', 'E', 'XW', 'X', 'x'); 
%-----------Multiplication--------------------%
XR=zeros(Np,P*L);
for k=1:P
    XR(:,(k-1)*L+1:k*L)=XD(:,k)*ones(1,L);
end
yc=ones(Np,1)*y';
XYM=XR.*yc;
XYM=conj(XYM');

% XD = conj(XD');                      %XDת�õĸ�����
% clear ('XF1', 'E', 'XW', 'X', 'x'); 
%-----------Multiplication--------------------%
% XM = zeros(P,Np^2);
% for k = 1:Np
%     for l = 1:Np
%         XM(:,(k-1)*Np+l) = (XD(:,k).*conj(XD(:,l)));
%     end
% end
% clear XD;

%------------Second FFT-----------------------%
XYF2 = fft(XYM);  
% clear XM;
XYF2 = fftshift(XYF2);
XYF2 = [XYF2(:,Np/2+1:Np^2) XYF2(:,1:Np^2/2)];
% length(XF2);
% XF2 = XF2(P/4:3*P/4,:);
M = abs(XYF2);
% clear XF2;  %Absolute value and complex magnitude

%%%%%%%%%%%%%%%%Ƶ�ʷֱ��ʺ�ѭ��Ƶ�ʷֱ���
alphao = (-1:1/N:1)*fs;     %about the variable N!!!
fo = (-0.5:1/Np:0.5)*fs;
Sxy = zeros(Np+1,2*N+1);     %about the variable N!!!
for k1 = 1:N
    for k2 = 1:Np
        alpha=(k1-1)/N+(k2-1)/Np-1;
        f=((k2-1)/Np-(k1-1)/N)/2;
        k=1+Np*(f+0.5);
        l=1+N*(alpha+1);
        Sxy(k,l)=M(k1,k2);
    end
end
Sxy = Sxy./max(max(Sxy)); % Normalizes the magnitudes of the values in output matrix (maximum = 1)
figure(1)
mesh(alphao,fo,Sxy)
title('SCD estimate using SSCA')
xlable('alpha');
ylable('f');
zlable('Sx');
figure(2)
contour(alphao,fo,Sxy);
xlabel('alpha(Hz)');
ylabel('f(Hz)');
figure(3)
plot(fo,Sxy(:,1+N*((alpha/fs+1))));
xlabel('f(Hz)');
ylabel('Sx(alpha)');