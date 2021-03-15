%% ������˼�룺��״�׹����㷨��SSCA��ʵ��ѭ�����ܶȹ���
%SSCA �㷨����һ�μ����������a=2fk-2fһ���������ź�ѭ���׵Ĺ��� 
clc
clear all
close all
%�������ò��� 
Fs=1280;
cf=160;
N=256;
PW=N/Fs;
M=64;                % �ɿ������� M=df/dalpha������Ԥ���趨������M��1
signal_type=3;
switch signal_type
    case 1; input_signal = sin(2*pi*(cf)*((0:(N-1))/Fs));
    case 2; input_signal = bpsk_generator(Fs,cf,N);
    case 3; input_signal = qpsk_generator(Fs,cf,N);
    case 4; input_signal = lfm_generator(Fs) ;  
%     case 4; input_signal = lfm_generator(Fs,cf,N);
%     case 5; input_signal = nlfm_generator(Fs*10^6,cf*10^6,PW*10^3);   
    otherwise;
end
% snr=10;
% x=awgn(input_signal,snr);

x=input_signal;

df=32;                  % ����Ƶ�ʷֱ���
N = (M*Fs)/df;
% N = M*N;
N = pow2 (nextpow2(N)); % ���ݹ۲ⳤ��

X = fft(x,N);           %���ź�Ƶ�׼�Ƶ�׹���
X = fftshift(X);
Xc = conj(X);                 

S = zeros (N,N);              % �洢ѭ�����ܶȹ���ֵ
f = zeros (N,N);              % �洢Ƶ��
alfa = zeros (N,N);           % �洢ѭ��Ƶ��
F = Fs/(2*N);                 % CSA��ԪƵ�ʼ��  F = Fs/(2*N);     
G = Fs/N;                     % CSA��Ԫѭ��Ƶ�ʼ�� -  G = Fs/N;  
m = -M/2+1:M/2;               % Ƶ��ƽ����������

for k = 1:N 
    %���ݹ۲ⳤ������Ƶ�ʺ�ѭ��Ƶ��
    k1 = 1:N;
    f(k,k1) = F*(k+k1-1) - Fs/2;         
    alfa(k,k1) = G*(k-k1 + N-1) - Fs;           
    for k1 = 1:N 
        B = max(1-k, 1-k1);          
        A = min (N-k, N-k1);         
        n = m((m<=A) & (m>=B));   %ͨ���ضϴ������޸�����������Χ������                                                 
        if isempty(n)
            S(k,k1) = 0;
        else
            p = k+n; q = k1+n;
            Y = X(p).*Xc(q);
            S(k,k1) = sum(Y);     %����ѭ�����ܶȵĹ���ֵ
        end
    end
end
Sx = abs(S./max(max(S)));     %��һ�����
figure(1)
mesh(alfa,f,Sx);
title('cyclic spectrum of BPSK signal')
xlabel('Frequency ')
ylabel('alpha ')
zlabel('magutide')
figure(2);
plot(alfa,Sx,'r'); 
xlabel('Cycle frequency (Hz)');  
ylabel('magutide')
figure(3);
plot(f,Sx);
xlabel('Frequency ')
ylabel('magutide')
% figure(1); 
% plot(SS)
% figure(2); 
% % f=zeros(N,N); 
% plot(alfa,Sx,'r'); 
% axis([-1000 1000 0 1]); 
% ylabel('Frequency (Hz)');xlabel('Cycle frequency (Hz)');  
% title('BPSKѭ������άͼ');
% figure(3)
% mesh(f,alfa,Sx);
% title('cyclic spectrum of BPSK signal')
% xlabel('Frequency ')
% ylabel('alpha ')
% zlabel('magutide')