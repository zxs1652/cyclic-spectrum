% function [Sx,alphao,fo] = autofam(x,fs,df,dalpha)%************************************************
% xΪ�����ź�����
% fsΪ������
% dfΪƵ�ʷֱ���
% dalphaΪѭ��Ƶ�ʷֱ���
% ע�⣺������df>>dalpha���ܵõ��ϸߵ�Ч����fs/dalphaΪ���в������㹦���׵���Ԫ����
%%fam�㷨��ssca(autossca)�㷨��,��Ч��һ��.������ѭ��ƽ������
%************************************************
% if nargin > 4 | nargin < 4
%     error('Wrong number of arguments.');
% end   

% fs=1280;
% fc=160;
% N=8192;

fs=1000;
fc=400;
N=5000;

PW=N/fs;
signal_type=2;
switch signal_type
    case 1; input_signal = sin(2*pi*(fc)*((0:(N-1))/fs));
    %case 2; input_signal = bpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 2; input_signal = bpsk_generator();
    case 3; input_signal = qpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 4; input_signal = lfm_generator(fs) ;
    case 5; input_signal = nlfm_generator(fs) ;
    otherwise;
end
x = awgn(input_signal,10);

%M=128;
M=16;
dalpha=fs/N;
df=dalpha*M;
%-----------Definition of Parameters-----------%
Np = pow2(nextpow2(fs/df));     %�����ŵ���
L = Np/4;                       %��ȡ���ӣ�ͬһ�������е����ڵ��ƫ�ƣ�L��ʾÿ�λ��������ݵ���
P = pow2(nextpow2(fs/dalpha/L));%�ŵ����������,P��ʾ��������
N = P*L;                        %�������ݵĵ���
%----------Input Channelization----------------%
if length(x) < N
    x(N) = 0;
elseif length(x) > N
    x = x(1:N);
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
%------------------Windowing------------------%
% a = hamming(Np);
a=kaiser(Np);
XW = diag(a)*X;   %ÿ�μӴ�,����Ƶ��й¶
% XW = X;         %���Ӵ�������Ӵ���Ч�����Ƚ�
%---------------��һ�θ���Ҷ�任--------------------%
XF1 = fft(XW);  
% clear XW;   
XF1 = fftshift(XF1);   
XF1 = [XF1(:,P/2+1:P) XF1(:,1:P/2)];  %�����еĹ����ǰѸ���Ҷ�任�Ľ���ϰ����°�齻��,������Ұ�黥��

%-------------�±�Ƶ------------------%
E = zeros(Np,P);

for k = -Np/2:Np/2-1
    for m = 0:P-1
        E(k+Np/2+1,m+1) = exp(-i*2*pi*k*m*L/Np);
    end
end
XD = XF1.*E;   
% clear XF1; 
XD = conj(XD');                      %XDת�õĸ�����
% clear ('XF1', 'E', 'XW', 'X', 'x'); 
%-----------�˻�����--------------------%
XM = zeros(P,Np^2);
for k = 1:Np
    for l = 1:Np
        XM(:,(k-1)*Np+l) = (XD(:,k).*conj(XD(:,l)));
    end
end
% clear XD;
%-----------�ڶ��θ���Ҷ�任-----------------------%
XF2 = fft(XM);  
% clear XM;
XF2 = fftshift(XF2);
XF2 = [XF2(:,Np^2/2+1:Np^2) XF2(:,1:Np^2/2)];
% length(XF2);
XF2 = XF2(P/4:3*P/4,:);
M = abs(XF2);
% clear XF2;  %Absolute value and complex magnitude

%%%%%%%%%%%%%%%%Ƶ�ʷֱ��ʺ�ѭ��Ƶ�ʷֱ���%%%%%%%%%%
alphao = (-1:1/N:1)*fs;     %about the variable N!!!
fo = (-0.5:1/Np:0.5)*fs;
Sx = zeros(Np+1,2*N+1);     %about the variable N!!!
for k1 = 1:P/2+1
    for k2 = 1:Np^2
        if rem(k2,Np) == 0
            l = Np/2-1;    
        else
            l = rem(k2,Np)-Np/2-1; 
        end
        k = ceil(k2/Np)-Np/2-1; 
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
            Sx(round(kk), round(ll)) = M(k1,k2); 
        end
    end
end
Sx = Sx./max(max(Sx)); % Normalizes the magnitudes of the values in output matrix (maximum = 1)
 
%figure 
%contour (alpha0, f0, Sx); grid;
%xlabel('Cycle frequency (Hz)'); ylabel('Frequency (Hz)');
%figure
%mesh (alpha0, f0, Sx); grid;
%xlabel('Cycle frequency (Hz)'); ylabel('Frequency (Hz)');
% title (['Time Smoothing SCD  ', filename, ', df = ', int2str(df),', N = ', int2str(N)]);
% %colorbar;
%a=max(max(Sx));     %��һ��
%Sx=Sx./a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ϊ�㷨��ͼ
figure(1);
mesh(alphao,fo,Sx);
title('SCD estimate using FAM')
xlabel('Cycle frequency(alpha)')
ylabel('frequency(f)')
zlabel('Sx')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
plot(alphao,Sx)
xlabel('Cycle frequency(alpha)')
ylabel('Sx')
figure(3);
plot(fo,Sx)
xlabel('frequency(f)')
ylabel('Sx')
% xlabel('Cycle frequency(alpha)')
% ylabel('frequency(f)')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3);
% plot(alphao,10*log10((Sx(Np/2+1,:))))%f=0�����δ��һ��,���߶���.
% xlabel('Cycle frequency(alpha)')
% title('f=0')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% figure(4);
% plot(fo,10*log10(Sx(:,N+1)))%alpha=0�����
% xlabel('frequency(f)')
% title('alpha=0')
% grid on;
% title('CMMB Signal');
% end
% figure(3)
% plot(alpha0,Sx)
% figure(4)
% plot(fo,Sx)