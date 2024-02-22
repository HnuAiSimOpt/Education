function k = T_Stiffness(E,NU,t,xe,ye,ID)
	if ID == 1
         D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
	elseif ID == 2
         D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2];
	end

	xy = zeros(3,2);
	for i = 1:3
		xy(i,1)=xe(i);
		xy(i,2)=ye(i);
	end

	k = zeros(6,6);
	[J B Ae]=T_B(xy);
	k = t*Ae*B'*D*B;
end