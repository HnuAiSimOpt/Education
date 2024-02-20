function [shapeq,dhdsq,dhdnq,dhdzq]=feisoq8(s,n,z)

si=[-1,1,1,-1,-1,1,1,-1];
ni=[-1,-1,1,1,-1,-1,1,1];
zi=[-1,-1,-1,-1,1,1,1,1];
h1=[1,-1,1,-1,1,-1,1,-1];
h2=[1,-1,-1,1,-1,1,1,-1];
h3=[1,1,-1,-1,-1,-1,1,1];
h4=[-1,1,-1,1,1,-1,1,-1];

for i=1:8
    
    shapeq(i)=0.125*(1+si(i)*s)*(1+ni(i)*n)*(1+zi(i)*z);
     dhdsq(i)=0.125*si(i)*(1+ni(i)*n)*(1+zi(i)*z);
     dhdnq(i)=0.125*ni(i)*(1+si(i)*s)*(1+zi(i)*z);
     dhdzq(i)=0.125*zi(i)*(1+si(i)*s)*(1+ni(i)*n)
     
end


end

