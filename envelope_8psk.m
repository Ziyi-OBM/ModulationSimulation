function [ I,Q ] = envelope_8psk( data )
%data is binary data (i.e 0 and 1) in a row vector
%Generate I and Q value for 8psk

I=zeros(1,length(data)/3);
Q=zeros(1,length(data)/3);
IQindex=1;
for n=1:3:length(data)
    b2=data(n);
    b1=data(n+1);
    b0=data(n+2);
     
    %Assign I
     
    if b1==1 && b0==0                  
        I(IQindex)=cos(pi/8);
        
    elseif b1==1 &&b0==1
        I(IQindex)=cos(3*pi/8);
        
    elseif b1==0 &&b0==1
        I(IQindex)=-cos(3*pi/8);
    elseif b1==0 &&b0==0
        I(IQindex)=-cos(pi/8);
    end
    
    %Assign Q
        if b2==1 && b0==0               
        Q(IQindex)=sin(pi/8);
        
    elseif b2==1 &&b0==1
        Q(IQindex)=sin(3*pi/8);
        
    elseif b2==0 &&b0==1
        Q(IQindex)=-sin(3*pi/8);
        
    elseif b2==0 &&b0==0
        Q(IQindex)=-sin(pi/8);
        
        end
        
        IQindex=IQindex+1;
   
end
        









end


