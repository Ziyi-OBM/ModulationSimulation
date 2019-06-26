function [ data ] = generatedata( n,an0,an1 )
%generatedata(number of data, value for 0, value for 1) 


bidata=round(rand(1,n));
data=zeros(1,n);

for i=1:n
    if bidata(i)==0
        data(i)=an0;
    else if bidata(i)==1
            data(i)=an1;
        end
    end
    
end

end

