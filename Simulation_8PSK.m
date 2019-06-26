clc;
close all;
clear all;

%define variables
D=100;  n=3*D;   %# of symbol / bits  
an0=0; %value assigned to 0
an1=1;  %value assigned to 1
kt=5;   %Each pulse last for 2*kt*tb seconds
R=7000;           %bitrate
bs=32  ;   %sample per bit
r=0.5;    %rolloff  factor
fc=10000;       %carrier frequency (Hz)
wc=2*pi*fc;
tb=1/R;   %bit period 0.0002
ts=tb/bs;       %sampleing time

fs=1/ts;        %sampleing frequency

noisestep=[0 0.01,0.04,0.25,0.49,1,2.25,4,];
sig=0;    % variance of noise

for noiseLvnum=1:length(noisestep)

sig=noisestep(noiseLvnum);

data=generatedata(n,an0,an1);   %generate data
% stem(data);
% title('Data points');



[ I,Q ] = envelope_8psk( data );     %Encode

I_length=length(I);
Q_length=length(I);

[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal
pause(0.5)

% figure
% plot(pt,h)
% title('plot the pluse')
% pause(0.5)

%generate the baseband OOK signal

t=-kt*tb:ts:tb*I_length+kt*tb;            %time vector of the signal

%generate 

I_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       I_signal(k)=I_signal(k)+I(m+1)*h(k-m*bs);       
       
    else
        
        I_signal(k)=I_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

Q_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       Q_signal(k)=Q_signal(k)+Q(m+1)*h(k-m*bs);       
       
    else
        
        Q_signal(k)=Q_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

clear k

% figure();
% subplot(2,1,1)
% plot(t,I_signal,0:tb:(I_length-1)*tb,I,'*');
% title('I signal using rRCRO')
% xlabel('time t')
% ylabel('I(t)')
% 
% subplot(2,1,2)
% plot(t,Q_signal,0:tb:(Q_length-1)*tb,Q,'*');
% title('Q signal using rRCRO')
% xlabel('time t')
% ylabel('Q(t)')
% 
% pause(0.5)





%modulation
% tmol=-kt*tb:ts:tb*n+kt*tb;
g1=cos(wc*t);
I_modulated=I_signal.*g1;

g2=sin(wc*t);
Q_modulated=Q_signal.*g2;

% figure();
% subplot(2,1,1)
% plot(t,I_modulated)
% title('Modulated signal I')
% xlabel('time t')
% ylabel('I(t)')
% 
% subplot(2,1,2)
% plot(t,Q_modulated)
% title('Modulated signal Q')
% xlabel('time t')
% ylabel('Q(t)')
% 
% pause(0.5)


TX_signal_nonoise=I_modulated-Q_modulated;


%add noise
% noise=random('normal',0,sqrt(sig),[1,length(tmol)]);
noise=randn(1,length(t))*sqrt(sig);

TX_signal=TX_signal_nonoise+noise;     %signal with noise

% figure();
% plot(t,TX_signal_nonoise);
% title('Modulated 8-PSK signal at the input of the channel')
% xlabel('time t')
% ylabel('s(t)')
% 
% pause(0.5)
% 
% figure();
% plot(t,TX_signal);
% title('Modulated 8-PSK signal at the output of the channel (with noise)')
% xlabel('time t')
% ylabel('s(t)')
% 
% pause(0.5)



%demodulation

I_demodulated=TX_signal.*cos(wc*t);             %demodulate I signal

t2=-2*kt*tb:ts:((length(I_demodulated)+length(h)-1)-1)*ts-2*kt*tb;    %correction for shift in time of -kt*tb
                                      
I_recovered=conv(I_demodulated,h);  %matched filter

Q_demodulated=TX_signal.*-sin(wc*t);             %demodulate Q signal

Q_recovered=conv(Q_demodulated,h);  %matched filter

%plotting demodulated signal
% 
% figure
% subplot(2,1,1);
% plot(t2,I_recovered,0:tb:(I_length-1)*tb,I*bs/2,'*')
% title('Recovered I signal')
% xlabel('t')
% ylabel('I(t)')
% 
% subplot(2,1,2);
% plot(t2,Q_recovered,0:tb:(I_length-1)*tb,Q*bs/2,'*')
% title('Recovered Q signal')
% xlabel('t')
% ylabel('Q(t)')
% 
% pause(0.5)


% calculate I signal at detection points
sample_time=0:tb:(D-1)*tb;
I_sampled=zeros(1,length(sample_time));
epsilon=tb/(10*bs);
for k1=0:D-1
for k2=1:length(t2)
    if t2(k2)>tb*k1-epsilon && t2(k2)<tb*k1+epsilon
    I_sampled(k1+1)=I_recovered(k2);
   
    end
end
end
clear k1
clear k2

% calculate Q signal at detection points
sample_time=0:tb:(D-1)*tb;
Q_sampled=zeros(1,length(sample_time));
epsilon=tb/(10*bs);
for k1=0:D-1
for k2=1:length(t2)
    if t2(k2)>tb*k1-epsilon && t2(k2)<tb*k1+epsilon
    Q_sampled(k1+1)=Q_recovered(k2);
   
    end
end
end
clear k1
clear k2



%detect error

data_detected=zeros(1,length(data));
IQindex=1;

for k=1:3:length(data)
    %     b2=data(n);
    %     b1=data(n+1);
    %     b0=data(n+2);
    
    if Q_sampled(IQindex)>0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=1;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=1;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=0;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=1;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=0;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=0;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=1;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=0;
        data_detected(k+2)=0;
        
    end
    
    IQindex=IQindex+1;
    
end

clear k

error_logical=data==data_detected;

correct_index=find(error_logical);

% correct_symbolindex=ceil(correct_index/3);
% correct_symbolindex=correct_symbolindex(1:3:length(correct_symbolindex));

error_index=find(~error_logical);

error_symbolindex=unique(ceil(error_index/3));


biterror_number=length(error_index);



%%%%plotting everything

if noiseLvnum==1

 stem(data);
title('Data points');

figure();
subplot(2,1,1)
plot(t,I_signal,0:tb:(I_length-1)*tb,I,'*');
title('I signal using rRCRO')
xlabel('time t')
ylabel('I(t)')

subplot(2,1,2)
plot(t,Q_signal,0:tb:(Q_length-1)*tb,Q,'*');
title('Q signal using rRCRO')
xlabel('time t')
ylabel('Q(t)')

pause(0.5)
figure();
subplot(2,1,1)
plot(t,I_modulated)
title('Modulated signal I')
xlabel('time t')
ylabel('I(t)')

subplot(2,1,2)
plot(t,Q_modulated)
title('Modulated signal Q')
xlabel('time t')
ylabel('Q(t)')

pause(0.5)
figure();
plot(t,TX_signal_nonoise);
title('Modulated 8-PSK signal at the input of the channel')
xlabel('time t')
ylabel('s(t)')

pause(0.5)



end


figure();
plot(t,TX_signal);
title('Modulated 8-PSK signal at the output of the channel (with noise)')
xlabel('time t')
ylabel('s(t)')

pause(0.5)

figure
subplot(2,1,1);
plot(t2,I_recovered,0:tb:(I_length-1)*tb,I*bs/2,'*')
title('Recovered I signal')
xlabel('t')
ylabel('I(t)')

subplot(2,1,2);
plot(t2,Q_recovered,0:tb:(I_length-1)*tb,Q*bs/2,'*')
title('Recovered Q signal')
xlabel('t')
ylabel('Q(t)')

pause(0.5)





figure
plot(I_sampled,Q_sampled,'*',I_sampled(error_symbolindex),Q_sampled(error_symbolindex),'r*')
grid on
title(['Received Constellation with noise of variance ',num2str(sig)])
xlabel('I')
ylabel('Q')


biterror_rate=biterror_number/n;



disp(['bit error rate is ',num2str(biterror_rate)])

if noiseLvnum<length(noisestep)

disp(['Press to simulate noise variance ',num2str(noisestep(noiseLvnum+1))])

end

pause

end


% disp(['bit error rate when variance is ',num2str(sig),' is ',num2str(errrate)])
 



%%%%%%%%%%%%%%%%%%%%%%find averaged power spectrum density%%%%%%%%%%%%




PSD_trials=100;  %repeat n time;
Etempt=zeros(1,PSD_trials);
psd_tempt=zeros(1,length(t)-(2*kt+1)*bs);

for k2=1:PSD_trials

data=generatedata(n,an0,an1);   %generate data

[ I,Q ] = envelope_8psk( data );     %Encode

I_length=length(I);
Q_length=length(I);

[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal

%generate the baseband OOK signal

t=-kt*tb:ts:tb*I_length+kt*tb;            %time vector of the signal

%generate 

I_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       I_signal(k)=I_signal(k)+I(m+1)*h(k-m*bs);       
       
    else
        
        I_signal(k)=I_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

Q_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 


for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       Q_signal(k)=Q_signal(k)+Q(m+1)*h(k-m*bs);       
       
    else
        
        Q_signal(k)=Q_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

clear k

%modulation
% tmol=-kt*tb:ts:tb*n+kt*tb;
g1=cos(wc*t);
I_modulated=I_signal.*g1;

g2=sin(wc*t);
Q_modulated=Q_signal.*g2;

TX_signal=I_modulated-Q_modulated;


TX_psd=zeros(1,length(TX_signal)-(2*kt+1)*bs);
t_psd=zeros(1,length(TX_signal)-(2*kt+1)*bs);
for n2=1:length(TX_psd)
    TX_psd(n2)=TX_signal(n2+kt*bs);
    t_psd(n2)=t(n2+kt*bs);
end

%find PSD (without noise)

NFFT=length(TX_psd);
sfft=fftshift(fft(TX_psd,NFFT))*ts;
f=2*(-(NFFT-1)/2:(NFFT-1)/2)./(NFFT-1)/(2*ts);

psd=sfft.*conj(sfft)/(tb*n);

psd_average=(psd+psd_tempt*(k2-1))/(k2);

psd_tempt=psd_average;

Etempt(k2)=trapz(t,TX_signal.^2);



end

%theoretical PSD

% p341= theoryPSD(f,fc,tb)   ;               
Eb=mean(Etempt)/n;

figure 
plot(t_psd,TX_psd)

figure;
plot(f,psd_average)
title(['Power spetral density of 8-PSK signal with bit rate ',num2str(R),' bit/s'])
xlabel('Frequency Hz')
ylabel('power')
grid on

pause(0.5)

figure;
plot(f,20.*log10(psd_average))
title(['Power spetral density in dB of 8-PSK signal with bit rate ',num2str(R),' bit/s'])
xlabel('Frequency Hz')
ylabel('power (dB)')
grid on







%%%%%%%%%%%%%%%%%%%%%%find bit error rate%%%%%%%%%%%%
pause

D=40;  n=3*D;   %redefine # of bits / symbol    

number_each_sig=50;  %repeat n time for each variance;
steps=30;            %steps of the BER plot
ENdb=linspace(-30,25,steps);            %define E/N to be measured in dB
EN=10.^(ENdb./10);                      %convert to normal
sigsq=Eb./(EN*2)*fs;                    %variance sigma^2
biterrplot=zeros(1,length(sigsq));
symerrplot=zeros(1,length(sigsq));
biterror_tempt=zeros(1,length(number_each_sig));
symerror_tempt=zeros(1,length(number_each_sig));

clear trials

for sigcount=1:length(sigsq);
    
for trials=1:number_each_sig
    

data=generatedata(n,an0,an1);   %generate data

[ I,Q ] = envelope_8psk( data );     %Encode

I_length=length(I);
Q_length=length(I);

[h,pt]=generateRRCRO(kt,tb,bs,r);    %create pulse used to generate signal

%generate the baseband OOK signal

t=-kt*tb:ts:tb*I_length+kt*tb;            %time vector of the signal

%generate 

I_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 

clear k

for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       I_signal(k)=I_signal(k)+I(m+1)*h(k-m*bs);       
       
    else
        
        I_signal(k)=I_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end

Q_signal=zeros(1,length(t));           %Signal waveform value correspond to the time 

clear k

for k=1:length(t)-1             %calculate signal
    
  for m=0:I_length-1                     %restrict n so that t-nTb is the domain of each pulse
    
    if (k-m*bs)>1 && (k-m*bs)<2*kt*bs        %only keep n's where h(t-nTb) is defined
        
       Q_signal(k)=Q_signal(k)+Q(m+1)*h(k-m*bs);       
       
    else
        
        Q_signal(k)=Q_signal(k)+0;            %if t-nTb is out side the domain, assign 0
    
    end
    
  end

end



%modulation
% tmol=-kt*tb:ts:tb*n+kt*tb;
g1=cos(wc*t);
I_modulated=I_signal.*g1;

g2=sin(wc*t);
Q_modulated=Q_signal.*g2;



TX_signal=I_modulated-Q_modulated;


%add noise
% noise=random('normal',0,sqrt(sig),[1,length(tmol)]);
noise=randn(1,length(t))*sqrt(sigsq(sigcount));

TX_signal=TX_signal+noise;     %signal with noise

%demodulation

I_demodulated=TX_signal.*cos(wc*t);             %demodulate I signal

t2=-2*kt*tb:ts:((length(I_demodulated)+length(h)-1)-1)*ts-2*kt*tb;    %correction for shift in time of -kt*tb
                                      
I_recovered=conv(I_demodulated,h);  %matched filter

Q_demodulated=TX_signal.*-sin(wc*t);             %demodulate Q signal

Q_recovered=conv(Q_demodulated,h);  %matched filter

%plotting demodulated signal

% calculate I signal at detection points
sample_time=0:tb:(D-1)*tb;
I_sampled=zeros(1,length(sample_time));
epsilon=tb/(10*bs);

clear k1
clear k2


for k1=0:D-1
for k2=1:length(t2)
    if t2(k2)>tb*k1-epsilon && t2(k2)<tb*k1+epsilon
    I_sampled(k1+1)=I_recovered(k2);
   
    end
end
end
clear k1
clear k2

% calculate Q signal at detection points
sample_time=0:tb:(D-1)*tb;
Q_sampled=zeros(1,length(sample_time));
epsilon=tb/(10*bs);
for k1=0:D-1
for k2=1:length(t2)
    if t2(k2)>tb*k1-epsilon && t2(k2)<tb*k1+epsilon
    Q_sampled(k1+1)=Q_recovered(k2);
   
    end
end
end




%detect error

data_detected=zeros(1,length(data));
IQindex=1;

clear k

for k=1:3:length(data)
    %     b2=data(n);
    %     b1=data(n+1);
    %     b0=data(n+2);
    
    if Q_sampled(IQindex)>0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=1;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=1;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=0;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=1;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)>0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=1;
        data_detected(k+1)=0;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))>abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=0;
        data_detected(k+2)=1;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)>0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=1;
        data_detected(k+2)=0;
    elseif Q_sampled(IQindex)<0 && I_sampled(IQindex)<0 && abs(Q_sampled(IQindex))<abs(I_sampled(IQindex))
        data_detected(k)=0;
        data_detected(k+1)=0;
        data_detected(k+2)=0;
        
    end
    
    IQindex=IQindex+1;
    
end



error_logical=data==data_detected;     %data_logical has 1 for correct and 0 for error in the cooresponding entry

correct_index=find(error_logical);

% correct_symbolindex=ceil(correct_index/3);
% correct_symbolindex=correct_symbolindex(1:3:length(correct_symbolindex));

error_index=find(~error_logical);

error_symbolindex=unique(ceil(error_index/3));

biterror_number=length(error_index);

biterror_rate=biterror_number/n;

symerror_number=length(error_symbolindex);

symerror_rate=symerror_number/D;
    
   
biterror_tempt(trials)=biterror_rate ;
symerror_tempt(trials)=symerror_rate ;

end

biterrplot(sigcount)=mean(biterror_tempt);
symerrplot(sigcount)=mean(symerror_tempt);

disp(['Step ',num2str(sigcount),' of ',num2str(length(sigsq))])




end



figure
plot(ENdb,biterrplot,ENdb,4/7*qfunc(sqrt((EN)*3*sin(pi/8)^2)));
title('Measured BER & Theoretical BER bound vs Eb/N0')
xlabel('E/N0')
ylabel('Pbe')
legend('Measured BER','Theoretical BER bound')

figure
plot(ENdb,symerrplot,ENdb,qfunc(sqrt((EN)*3*sin(pi/8)^2)));
title('Measured & Theoretical symbol error rate vs Eb/N0')
xlabel('E/N0')
ylabel('Pde')
legend('symbol error rate','Theoretical symbol error rate bound')





