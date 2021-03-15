%% ������˼�룺ʱ��ƽ��FFT�ۼ��㷨(FAM)����ѭ�����ܶ�
fs=1280;
fc=160;
N=8192;
PW=N/fs;

signal_type=1;
switch signal_type
    case 1; input_signal = sin(2*pi*(fc)*((0:(N-1))/fs));
    case 2; input_signal = bpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 3; input_signal = qpsk_generator(fs*10^6,fc*10^6,PW*10^3);
    case 4; input_signal = lfm_generator(fs*10^6,fc*10^6,PW*10^3);
    otherwise;
end

sig = awgn(input_signal,0);

M=128;
dalpha=fs/N;
df=dalpha*M;

% ��������
Np = pow2(nextpow2(fs/df));         % ��FFT����N'                                    
L = Np/4;                           % ���ڶ�ʱFFT����ص����ӣ�L<N'/4                                    
P = pow2(nextpow2(fs/dalpha/L));    % �ڶ���FFT�任����
N = P*L;                            % �ܹ۲ⳤ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����źŴ���
if length(sig)<N          
    sig(N) = 0; 
elseif length (sig)>N
    sig = sig(1:N);
end

NN = (P-1)*L+Np;
signal = sig;
signal(NN) = 0;
signal = signal(:);
X = zeros(Np,P);
for k = 0:P-1
    X(:,k+1) = signal(k*L+1:k*L+Np);    % �źŷֶδ���
end

a = hamming (Np);
XW = diag(a)*X;         %ѡȡhamming����ø��õĸ��ǽ��

XF1 = fft(XW); 
XF1 = fftshift(XF1);
XF1 = [XF1(:,P/2+1:P) XF1(:,1:P/2)];        %�źŵ�һ��FFT�任

%% �±�Ƶ��ø��������
E = zeros(Np,P);
for k = -Np/2:Np/2-1
    for m = 0:P-1
        E(k+Np/2+1, m+1) = exp(-i*2*pi*k*m*L/Np);       
    end
end

XD = XF1.*E;        %�����
XD = conj(XD');     %��ȡ������

clear ('XF1', 'E', 'XW', 'X', 'x','I', 'Q'); 

%% ��������г˻�
XM = zeros(P,Np^2);
for k = 1:Np
    for c = 1:Np
        XM(:,(k-1)*Np+c) = (XD(:,k).*conj(XD(:,c)));
    end
end

%%  �ڶ���P��FFT����

XF2 = fft(XM);
XF2 = fftshift(XF2);
XF2 = [XF2(:,Np^2/2+1:Np^2) XF2(:,1:Np^2/2)];
XF2 = XF2 (P/4:3*P/4,:);
MM = abs(XF2);

alpha0 = -fs:fs/N:fs;   
f0 = -fs/2:fs/Np:fs/2;  
Sx = zeros(Np+1, 2*N+1);    

for k1 = 1:P/2+1
    for k2 = 1:Np^2
        if rem(k2,Np) == 0
            c = Np/2 - 1;
        else
            c = rem(k2,Np) - Np/2 - 1;
        end
        k = ceil(k2/Np)-Np/2-1;
        p = k1-P/4-1;
        alpha = (k-c)/Np+(p-1)/N;
        f = (k+c)/2/Np;
        if alpha<-1 | alpha>1
            k2 = k2+1;
        elseif f<-.5 | f>.5
            k2 = k2+1;
        else
            kk = 1+Np*(f + .5);
            ll = 1+N*(alpha + 1);
            Sx(round(kk), round(ll)) = MM(k1,k2);
        end
    end
end
clear ('alpha', 'XM', 'XF2', 'MM', 'f'); 

Sx = Sx./max(max(Sx)); % ��һ�����
 
 figure 
 mesh (alpha0, f0, Sx); grid;
 xlabel('Cycle frequency (Hz)'); ylabel('Frequency (Hz)');
figure
plot(alpha0,Sx(65,:))
% mesh(f,alpha0,Sx)