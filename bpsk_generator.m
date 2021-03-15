function psk=bpsk_generator()

i=10;
j=5000;
fc=400;
% fm=i/5;
% B=2*fm
t=linspace(0,5,j);
a=round(rand(1,i));
st1=t;
for n=1:10;
    if a(n)<1
        for m=j/i*(n-1)+1:j/i*n;
            st1(m)=0;
        end
    else
        for m=j/i*(n-1)+1:j/i*n;
            st1(m)=1;
        end
    end
end
st2=t;
for k=1:j;
    if st1(k)>=1
        st2(k)=0;
    else
        st2(k)=1;
    end
end

st3=st1-st2;
% figure()
% plot(t,st3)

s1=cos(2*pi*fc*t);
psk=st3.*s1;
% figure()
% plot(t,psk)
% axis([0,1,-1,2])

end
% x=psk;
% N0=5000
% nn=0:length(x)-1
% for k=0:N0-1
%     xk(k+1)=sum(x.*exp(-j*2*pi*nn*k/N0));
% end
% figure()
% plot(xk)
% xk
