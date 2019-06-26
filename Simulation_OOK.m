clc;
close all;
clear all;

%define variables
n=32;   %# of bits
%n=100; %# of bits for error rate calculation
an0=0; %value assigned to 0
an1=1;  %value assigned to 1
kt=5;   %Each pulse last for 2*kt*tb seconds
tb=0.001;   %bit period
%bs=16;  %samples per bit
% ts=tb/bs;  %sampleing time
bs0=16  ;   %sample per carrier period
r=0.5;    %rolloff  factor
fc=10000;    %carrier frequency (Hz)
wc=2*pi*fc;
ts=1/(fc*bs0);       %sampleing time
bs=tb/ts;     %samples per bit
R=1/tb;  %bitrate
fs=bs/tb;        %sampleing frequency
sig=0.2;    % variance of noise


data=generatedata(n,an0,an1);   %generate data
stem(data);
title('Data points');

[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal

figure
plot(pt,h)
title('plot the pluse')

%generate the baseband OOK signal

t=-kt*tb:ts:tb*n+kt*tb;            %time vector of the signal

s=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:n-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       s(k)=s(k)+data(m+1)*h(k-m*bs);       
       
    else
        
        s(k)=s(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

figure();
plot(t,s,0:tb:(n-1)*tb,data,'*');
title(['RRCRO pulse baseband signal with r = ',num2str(r)])
xlabel('time t')
ylabel('signal amplitude s(t)')



%AM modulation
tmol=-kt*tb:ts:tb*n+kt*tb;
g=cos(wc*tmol);
modu=s.*g;

figure();
plot(tmol,modu);
title('Modulated signal ')
xlabel('time t')
ylabel('signal amplitude s(t)')


%add noise
% noise=random('normal',0,sqrt(sig),[1,length(tmol)]);
noise=randn(1,length(tmol))*sqrt(sig);

modun=modu+noise;     %signal with noise
figure();
plot(tmol,modun);
title('Modulated signal with noise')
xlabel('time t')
ylabel('signal amplitude s(t)')

%demodulation

Tx_mod=modun.*g;
NFFT=length(Tx_mod);
sfft=fftshift(fft(Tx_mod,NFFT))*ts;
f=2*(-(NFFT-1)/2:(NFFT-1)/2)./(NFFT-1)/(2*ts);


t2=-2*kt*tb:ts:((length(Tx_mod)+length(h)-1)-1)*ts-2*kt*tb;    %correction for shift in time of -kt*tb
                                      

RX=conv(Tx_mod,h);  %matched filter


%calculate signal at detection points
node=0:tb:(n-1)*tb;
RXnode=zeros(1,length(node));
epsilon=tb/(10*bs);
for kk=0:n-1
for k2=1:length(t2)
    if t2(k2)>tb*kk-epsilon && t2(k2)<tb*kk+epsilon
    RXnode(kk+1)=RX(k2);
   
    end
end
end
figure()
plot(t2,RX,0:tb:(n-1)*tb,data*bs/2,'r*',0:tb:(n-1)*tb,RXnode,'bo',t2,bs/4*ones(1,length(t2)),'b:');               
title(['RX baseband signal recovered , added by noise of variance ',num2str(sig),' in the channel'])
xlabel('time t')
ylabel('signal amplitude s(t)')

%error rate
err=0;

for w=1:n
    if data(w)==an1 && RXnode(w)<bs/4
        err=err+1;
        
    else if data(w)==an0 && RXnode(w)>bs/4
        err=err+1;
        
        end
    end

end



errrate=err/n;

disp(['bit error rate when variance is ',num2str(sig),' is ',num2str(errrate)])
 
pause



%%%%%%%%%%%%%%%%%%%%%%find averaged power spectrum density%%%%%%%%%%%%




pp=250;  %repeat n time;
Etempt=zeros(1,pp);
psp=zeros(1,length(tmol)-(2*kt+1)*bs+1);

for k2=1:pp


data=generatedata(n,an0,an1);   %generate data
% stem(data);
% title('Data points');


[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal




%generate the baseband OOK signal

t=-kt*tb:ts:tb*n+kt*tb;            %time vector of the signal

s=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:n-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       s(k)=s(k)+data(m+1)*h(k-m*bs); 
       
       
    else
        
        s(k)=s(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end





%AM modulation
tmol=-kt*tb:ts:tb*n+kt*tb;
g=cos(wc*tmol);

modu=s.*g;


sf=zeros(1,length(modu)-(2*kt+1)*bs+1);
t0=zeros(1,length(tmol)-(2*kt+1)*bs+1);
for n2=1:length(sf)
    sf(n2)=modu(n2+kt*bs);
    t0(n2)=tmol(n2+kt*bs);
end
modu=sf;
tmol=t0;

%find PSD (without noise)

NFFT=length(modu);
sfft=fftshift(fft(modu,NFFT))*ts;
f=2*(-(NFFT-1)/2:(NFFT-1)/2)./(NFFT-1)/(2*ts);

ps=sfft.*conj(sfft)/(tb*n);

psa=(ps+psp*(k2-1))/(k2);

psp=psa;

Etempt(k2)=trapz(tmol,modu.^2);



end

%theoretical PSD

p341= theoryPSD(f,fc,tb)   ;               
Eb=mean(Etempt)/n;



figure;
plot(f,psp,f,p341)
title(['Averaged power spectrum density over ',num2str(pp),' loops in dB'])
xlabel('Frequency Hz')
ylabel('power')
figure;
plot(f,20.*log10(psa),f,20.*log10(p341))
title(['Averaged power spectrum density over ',num2str(pp),' loops'])
xlabel('Frequency Hz')
ylabel('power (dB)')






%%%%%%%%%%%%%%%%%%%%%%find bit error rate%%%%%%%%%%%%
pause

nn=50;  %repeat n time for each variance;
ENdb=linspace(-20,25,30);
EN=10.^(ENdb./10);
sigsq=Eb./(EN*2)*fs;
errplot=zeros(1,length(sigsq));
etempt=zeros(1,length(nn));



for sigcount=1:length(sigsq);
    
for rep=1:nn
    
data=generatedata(n,an0,an1);   %generate data

[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal

%generate the baseband OOK signal

t=-kt*tb:ts:tb*n+kt*tb;            %time vector of the signal

s=zeros(1,length(t));           %Signal waveform value correspond to the time 

for k=1:length(t)-1             %calculate signal
    
  for m=0:n-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       s(k)=s(k)+data(m+1)*h(k-m*bs);       
       
    else
        
        s(k)=s(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

%AM modulation
% tmol=-kt*tb:ts:tb*n+kt*tb;
g=cos(wc*t);
modu=s.*g;

%add noise

noise=randn(1,length(t))*sqrt(sigsq(sigcount));

modun=modu+noise;     %signal with noise


%demodulation

Tx_mod=modun.*g;

t2=-2*kt*tb:ts:((length(Tx_mod)+length(h)-1)-1)*ts-2*kt*tb;    %correction for shift in time of -kt*tb
                                     
RX=conv(Tx_mod,h);                  %matched filter RX is the received signal

node=0:tb:(n-1)*tb;
RXnode=zeros(1,length(node));
epsilon=tb/(10*bs);
for kk=0:n-1
    for k2=1:length(t2)
        if t2(k2)>tb*kk-epsilon && t2(k2)<tb*kk+epsilon
            RXnode(kk+1)=RX(k2);
            
        end
    end
end


%error rate
err=0;

for w=1:n
    if data(w)==an1 && RXnode(w)<bs/4
        err=err+1;
    else if data(w)==an0 && RXnode(w)>bs/4
        err=err+1;
        end
    end

end
errrate=err/n;
etempt(rep)=errrate ;
end

errplot(sigcount)=mean(etempt);

disp(['Step ',num2str(sigcount),' of ',num2str(length(sigsq))])

end



figure
plot(ENdb,(errplot),ENdb,(qfunc(sqrt(EN))));
title('Pe vs SNR')
xlabel('E/N0')
ylabel('Pe')


