function ps = theoryPSD( f,fc,tb)
% [ pulsevlaue,timevector ] = generateRCRO( 1/2pulselength/tb,bitperiod,sampleperbit,roll of factor r )




R=1/tb;     %bitrate



ps=zeros(1,length(f));


for k=1:length(f)
    
   if abs(f(k)-fc)<R/4 
       ps(k)=tb/16*(1+dirac(f(k)-fc));
       
       elseif   abs(f(k)+fc)<R/4
           ps(k)=tb/16*(1+dirac(f(k)+fc));
      
        elseif abs(f(k)-fc)<(3*R/4) && abs(f(k)-fc)>=(R/4) 
           ps(k)=tb/32*(1+cos((2*pi*(abs(f(k)-fc)-R/4))/R));
           
        elseif abs(f(k)+fc)<(3*R/4) && abs(f(k)+fc)>=(R/4) 
           ps(k)=tb/32*(1+cos((2*pi*(abs(f(k)+fc)-R/4))/R));
           
       else
           ps(k)=0;     
   end
   


% plot(t,h)

end