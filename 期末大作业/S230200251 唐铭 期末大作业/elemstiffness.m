function ke=elemstiffness(E,nu,t,nodeloc)
%% node loc
x1=nodeloc(1,1);                
y1=nodeloc(1,2);
x2=nodeloc(2,1);                
y2=nodeloc(2,2);
x3=nodeloc(3,1);                
y3=nodeloc(3,2);
x4=nodeloc(4,1);                
y4=nodeloc(4,2);

x5=(x1+x2)/2;
y5=(y1+y2)/2;
x6=(x2+x3)/2;
y6=(y2+y3)/2;
x7=(x3+x4)/2;
y7=(y3+y4)/2;
x8=(x4+x1)/2;
y8=(y4+y1)/2;

%% shape function 
syms a b

N1=(1-a)*(1-b)*(-a-b-1)/4;
N2=(1+a)*(1-b)*(a-b-1)/4;
N3=(1+a)*(1+b)*(a+b-1)/4;
N4=(1-a)*(1+b)*(-a+b-1)/4;
N5=(1-b)*(1+a)*(1-a)/2;
N6=(1+a)*(1+b)*(1-b)/2;
N7=(1+a)*(1+b)*(1-a)/2;
N8=(1-a)*(1+b)*(1-b)/2;

x=N1*x1+N2*x2+N3*x3+N4*x4+N5*x5+N6*x6+N7*x7+N8*x8;
y=N1*y1+N2*y2+N3*y3+N4*y4+N5*y5+N6*y6+N7*y7+N8*y8;

dx1=diff(x,a);
dx2=diff(x,b);
dy1=diff(y,a);
dy2=diff(y,b);

N1a=diff(N1,a);
N1b=diff(N1,b);
N2a=diff(N2,a);
N2b=diff(N2,b);
N3a=diff(N3,a);
N3b=diff(N3,b);
N4a=diff(N4,a);
N4b=diff(N4,b);
N5a=diff(N5,a);
N5b=diff(N5,b);
N6a=diff(N6,a);
N6b=diff(N6,b);
N7a=diff(N7,a);
N7b=diff(N7,b);
N8a=diff(N8,a);
N8b=diff(N8,b);

%% B matrix
B11=dy2*N1a-dy1*N1b;
B13=dy2*N2a-dy1*N2b;
B15=dy2*N3a-dy1*N3b;
B17=dy2*N4a-dy1*N4b;
B19=dy2*N5a-dy1*N5b;
B111=dy2*N6a-dy1*N6b;
B113=dy2*N7a-dy1*N7b;
B115=dy2*N8a-dy1*N8b;

B12=0;
B14=0;
B16=0;
B18=0;
B110=0;
B112=0;
B114=0;
B116=0;

B22=dx1*N1b-dx2*N1a;
B24=dx1*N2b-dx2*N2a;
B26=dx1*N3b-dx2*N3a;
B28=dx1*N4b-dx2*N4a;
B210=dx1*N5b-dx2*N5a;
B212=dx1*N6b-dx2*N6a;
B214=dx1*N7b-dx2*N7a;
B216=dx1*N8b-dx2*N8a;

B21=0;
B23=0;
B25=0;
B27=0;
B29=0;
B211=0;
B213=0;
B215=0;


B31=dx1*N1b-dx2*N1a;
B32=dy2*N1a-dy1*N1b;
B33=dx1*N2b-dx2*N2a;
B34=dy2*N2a-dy1*N2b;
B35=dx1*N3b-dx2*N3a;
B36=dy2*N3a-dy1*N3b;
B37=dx1*N4b-dx2*N4a;
B38=dy2*N4a-dy1*N4b;
B39=dx1*N5b-dx2*N5a;
B310=dy2*N5a-dy1*N5b;
B311=dx1*N6b-dx2*N6a;
B312=dy2*N6a-dy1*N6b;
B313=dx1*N7b-dx2*N7a;
B314=dy2*N7a-dy1*N7b;
B315=dx1*N8b-dx2*N8a;
B316=dy2*N8a-dy1*N8b;


B=[B11 B12 B13 B14 B15 B16 B17 B18 B19 B110 B111 B112 B113 B114 B115 B116
   B21 B22 B23 B24 B25 B26 B27 B28 B29 B210 B211 B212 B213 B214 B215 B216 
   B31 B32 B33 B34 B35 B36 B37 B38 B39 B310 B311 B312 B313 B314 B315 B316];



A=1;
    
D=E/(1-nu^2)*[1 nu 0;
              nu 1 0;
              0 0 (1-nu)/2];
J=dx1*dy2-dy1*dx2;
B=simplify(B);
J=simplify(J);
k=A*B'*D*B/J; 

k1=int(int(k,b,-1,1),a,-1,1);
ke=t*k1;
ke=double(ke);


end