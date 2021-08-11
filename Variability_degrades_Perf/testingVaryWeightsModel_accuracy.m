function [p_ra]=testingVaryWeightsModel_accuracy(C,Wop,classAction1,numTrials,numtrainingSamples,Y)

              
numm=numTrials-numtrainingSamples;
p_ra=0;

i=1; 

odors= (size(Y,2))*numm;

for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(Y,2))

                        z1=Wop(:,1)'*Y(:,odour,trials);
                        z2=Wop(:,2)'*Y(:,odour,trials);

                       
                        
                        if (~ isempty(find(classAction1==odour)) )
                            
                            
                            pr_action1=  1/(1+exp(C*(z2-z1)));
                            
                            p_ra=p_ra+(pr_action1/odors);

                        end

                        if ( isempty((find(classAction1==odour))) )
                            
                            
                            pr_action2=  1/(exp(C*(z1-z2))+1);
                            
            
                            p_ra=p_ra+(pr_action2/odors);

                        end

                        
                    end
end



end