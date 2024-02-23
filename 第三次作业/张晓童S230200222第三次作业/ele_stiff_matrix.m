function ke=ele_stiff_matrix(ele,node,i,A,D,t)
      
      b1=node(ele(i,2),2)-node(ele(i,3),2); %b1=y1-y3
      b2=node(ele(i,3),2)-node(ele(i,1),2);
      b3=node(ele(i,1),2)-node(ele(i,2),2);
      c1=node(ele(i,3),1)-node(ele(i,2),1);%c1=x3-x2
      c2=node(ele(i,1),1)-node(ele(i,3),1);
      c3=node(ele(i,2),1)-node(ele(i,1),1);
      
      B(:,:,i)=1/(2*A)*[b1, 0 , b2, 0 , b3, 0 ;
                        0 , c1, 0 , c2, 0 , c3;
                        c1, b1, c2, b2, c3, b3];%单元B矩阵

      ke=B(:,:,i)'*D*B(:,:,i)*t*A; %单元刚度矩阵 