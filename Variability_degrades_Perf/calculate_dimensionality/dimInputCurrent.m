function dim_h= dimInputCurrent(C)


M=size(C,1);
Cii= C(find(eye(M)));

Cii_avg=mean(Cii(:));
Cii_varr= var(Cii(:));

Cij=C(find(~eye(M)));

Cij_avg=mean(Cij(:));
Cij_varr=var(Cij(:));

dim_h= (M*(Cii_avg^2))/( (Cii_avg^2)+ Cii_varr+ ((M-1)* ((Cij_avg^2)+(Cij_varr)) )         );

end