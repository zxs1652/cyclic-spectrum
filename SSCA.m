%% 本程序思想：带状谱估计算法（SSCA）实现循环谱密度估计
%SSCA 算法可以一次计算出，沿着a=2fk-2f一个条带上信号循环谱的估计 
clc
clear all
close all
%数据设置部分 
Fs=1280;
cf=160;
N=256;
PW=N/Fs;
M=64;                % 可靠性条件 M=df/dalpha，这里预先设定，满足M》1
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

df=32;                  % 定义频率分辨率
N = (M*Fs)/df;
% N = M*N;
N = pow2 (nextpow2(N)); % 数据观测长度

X = fft(x,N);           %求信号频谱及频谱共轭
X = fftshift(X);
Xc = conj(X);                 

S = zeros (N,N);              % 存储循环谱密度估计值
f = zeros (N,N);              % 存储频率
alfa = zeros (N,N);           % 存储循环频率
F = Fs/(2*N);                 % CSA单元频率间隔  F = Fs/(2*N);     
G = Fs/N;                     % CSA单元循环频率间隔 -  G = Fs/N;  
m = -M/2+1:M/2;               % 频率平滑窗口索引

for k = 1:N 
    %根据观测长度设置频率和循环频率
    k1 = 1:N;
    f(k,k1) = F*(k+k1-1) - Fs/2;         
    alfa(k,k1) = G*(k-k1 + N-1) - Fs;           
    for k1 = 1:N 
        B = max(1-k, 1-k1);          
        A = min (N-k, N-k1);         
        n = m((m<=A) & (m>=B));   %通过截断窗口来修复索引超出范围的问题                                                 
        if isempty(n)
            S(k,k1) = 0;
        else
            p = k+n; q = k1+n;
            Y = X(p).*Xc(q);
            S(k,k1) = sum(Y);     %计算循环谱密度的估计值
        end
    end
end
Sx = abs(S./max(max(S)));     %归一化输出
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
% title('BPSK循环谱三维图');
% figure(3)
% mesh(f,alfa,Sx);
% title('cyclic spectrum of BPSK signal')
% xlabel('Frequency ')
% ylabel('alpha ')
% zlabel('magutide')