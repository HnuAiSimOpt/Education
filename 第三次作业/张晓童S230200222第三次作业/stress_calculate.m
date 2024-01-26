function sgm=stress_calculate(ele,node,i,A,D,~,a)

      b1=node(ele(i,2),2)-node(ele(i,3),2);
      b2=node(ele(i,3),2)-node(ele(i,1),2);
      b3=node(ele(i,1),2)-node(ele(i,2),2);
      c1=-node(ele(i,2),1)+node(ele(i,3),1);
      c2=-node(ele(i,3),1)+node(ele(i,1),1);
      c3=-node(ele(i,1),1)+node(ele(i,2),1);
      
      B(:,:,i)=1/(2*A)*[b1,0,b2,0 ,b3,0;0,c1,0,c2,0,c3;c1,b1,c2,b2,c3,b3];
      ae(1,i)=a(2*ele(i,1)-1);
      ae(2,i)=a(2*ele(i,1));
      ae(3,i)=a(2*ele(i,2)-1);
      ae(4,i)=a(2*ele(i,2));
      ae(5,i)=a(2*ele(i,3)-1);
      ae(6,i)=a(2*ele(i,3));
      yps(:,i)=B(:,:,i)*ae(:,i);%单元应变
      sgm=D*yps(:,i);%单元应力
  end
  