%% STRESS SOLUTION
function stress=stresssolu(E,nu,nodeloc,elemdisp)
x1=nodeloc(1,1);                
y1=nodeloc(1,2);
x2=nodeloc(2,1);                
y2=nodeloc(2,2);
x3=nodeloc(3,1);                
y3=nodeloc(3,2);

A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;

b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;

B=0.5/A*[b1 0 b2 0 b3 0;
         0 c1 0 c2 0 c3;
         c1 b1 c2 b2 c3 b3];

D=E/(1-nu^2)*[1 nu 0;
              nu 1 0;
              0 0 (1-nu)/2];
     
%D=E*(1-nu)/(1-nu^2)/(1+nu)*[1 nu/(1-nu) 0;
%                            nu/(1-nu) 1 0;
%                            0 0 (1-2*nu)/2/(1-nu)];

stress=D*B*elemdisp;
          
end