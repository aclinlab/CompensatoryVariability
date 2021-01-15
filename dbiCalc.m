function PW_dbi=dbiCalc(Y)

no_ods= 100;
PW_dbi=zeros(no_ods);
numTrials=30;
n=2000;
Y=rescale(Y);

%normalize the vector responses to have a unity magnitude
Y=Y./ repmat(vecnorm(Y),n,1);

Y= reshape(Y,n,no_ods,numTrials);

for i=1:no_ods
    
    S1Centroid=mean(Y(:,i,:),3);
    
    Yreshaped= permute(Y(:,i,:),[1,3,2]);
    
    S1Var=sqrt(mean(  vecnorm((Yreshaped-S1Centroid ))    ));
           
    
    for j=1:no_ods
        
        S2Centroid=mean(Y(:,j,:),3);
    
        Yjreshaped= permute(Y(:,j,:),[1,3,2]);
    
        S2Var=sqrt(mean(  vecnorm((Yjreshaped-S2Centroid ))    ));
           
        PW_dbi(i,j)= (S1Var+S2Var)/(norm(S1Centroid-S2Centroid));
                
        
    end
    
    
    
end








end