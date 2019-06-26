function [ h,t ] = generateRRCRO( kt,tb,bs,r )
% [ pulsevlaue,timevector ] = generateRCRO( 1/2pulselength/tb,bitperiod,sampleperbit,roll of factor r )

%%%%test
% kt=5;
% tb=1;
% bs=16;
% r=0.8;
%%%%%%%%

ts=tb/bs;
R=1/tb;     %bitrate

t=-kt*tb:ts:kt*tb;

h=zeros(1,length(t));

epsilon=tb/(10*bs);

for k=1:length(t)
    
   if t(k)==0
       h(k)=1-r+4*r/pi;
      
   elseif abs(t(k))<(tb/(4*r)+epsilon) &&  abs(t(k))>(tb/(4*r)-epsilon)
           h(k)=r/sqrt(2)*((1+2/pi)*sin(pi/(4*r))+(1-2/pi)*cos(pi/(4*r)));
           
       else
           h(k)=(sin(pi*R*t(k)*(1-r))+4*R*r*t(k)*cos(pi*R*t(k)*(1+r)))/(pi*R*t(k)*(1-(4*R*r*t(k))^2));           
      
   end
    
end


% plot(t,h)

end



