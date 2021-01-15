function PW_AngDist=AngDist_Pw_Calc(Y)
% y is reshaped: 2000 x odorsx num trials
% y is rescaled from 0 to 1
n=2000;
odors=100;
numTrials=30;

Y= reshape(Y,n,odors,numTrials);
Y=rescale(Y);

no_ods= size(Y,2);
PW_AngDist=zeros(no_ods).*inf;

for i=1:no_ods
    
    S1Centroid=mean(Y(:,i,:),3);
    
    Yreshaped= permute(Y(:,i,:),[1,3,2]);               
    
    for j=1:no_ods
        
        if(i~=j)
        
            S2Centroid=mean(Y(:,j,:),3);

            Yjreshaped= permute(Y(:,j,:),[1,3,2]);

            PW_AngDist(i,j)= (acos(( (S1Centroid'/norm(S1Centroid'))*(S2Centroid/norm(S2Centroid))) )/ (0.5*pi) ) ;
        end
        
    end
    
    
    
end

end



