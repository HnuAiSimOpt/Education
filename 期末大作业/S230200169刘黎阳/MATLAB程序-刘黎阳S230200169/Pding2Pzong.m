function P9=Pding2Pzong(p4,JMV)

for i=1:length(JMV(:,1))

    pp=[p4(JMV(i,1)),p4(JMV(i,3)),p4(JMV(i,9)),p4(JMV(i,7))];
    P9(JMV(i,1))=p4(JMV(i,1));
    kesi=0;ita=-1;
    Fyp= [1/4*(1-kesi)*(1-ita)
          1/4*(1+kesi)*(1-ita)
          1/4*(1+kesi)*(1+ita)
          1/4*(1-kesi)*(1+ita)];
    P9(JMV(i,2))=pp*Fyp;
    P9(JMV(i,3))=p4(JMV(i,3));
    kesi=-1;ita=0;
    Fyp= [1/4*(1-kesi)*(1-ita)
          1/4*(1+kesi)*(1-ita)
          1/4*(1+kesi)*(1+ita)
          1/4*(1-kesi)*(1+ita)];   
    P9(JMV(i,4))=pp*Fyp;
    kesi=0;ita=0;
    Fyp= [1/4*(1-kesi)*(1-ita)
          1/4*(1+kesi)*(1-ita)
          1/4*(1+kesi)*(1+ita)
          1/4*(1-kesi)*(1+ita)]; 

    P9(JMV(i,5))=pp*Fyp;
    kesi=1;ita=0;
    Fyp= [1/4*(1-kesi)*(1-ita)
          1/4*(1+kesi)*(1-ita)
          1/4*(1+kesi)*(1+ita)
          1/4*(1-kesi)*(1+ita)];

    P9(JMV(i,6))=pp*Fyp;
    P9(JMV(i,7))=p4(JMV(i,7));

    kesi=0;ita=1;

    Fyp= [1/4*(1-kesi)*(1-ita)
          1/4*(1+kesi)*(1-ita)
          1/4*(1+kesi)*(1+ita)
          1/4*(1-kesi)*(1+ita)];
    P9(JMV(i,8))=pp*Fyp;
    P9(JMV(i,9))=p4(JMV(i,9));
end