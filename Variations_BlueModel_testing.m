function [P_raEq2_noxjk,p_raEq2_wHy]=Variations_BlueModel_testing(C,WopEq2_noxjk,WopEq2_wHy,...
    classAction1,numTrials,numtrainingSamples,Y_comp2_noxjk,Y_comp2_wHy)

              
numm=numTrials-numtrainingSamples;
P_raEq2_noxjk=0;
p_raEq2_wHy=0;
odors= (size(Y_comp2_noxjk,2))*numm;

for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(Y_comp2_noxjk,2))

                    
                        z1Eq2_noxjk=WopEq2_noxjk(:,1)'*Y_comp2_noxjk(:,odour,trials);
                        z2Eq2_noxjk=WopEq2_noxjk(:,2)'*Y_comp2_noxjk(:,odour,trials);
                        
                        z1Eq2_wHy=WopEq2_wHy(:,1)'*Y_comp2_wHy(:,odour,trials);
                        z2Eq2_wHy=WopEq2_wHy(:,2)'*Y_comp2_wHy(:,odour,trials);

                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1Eq2_noxjk= 1/(1+exp(C*(z2Eq2_noxjk-z1Eq2_noxjk)));
                            P_raEq2_noxjk=P_raEq2_noxjk+(pr_action1Eq2_noxjk/odors);
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2Eq2_noxjk=  1/(exp(C*(z1Eq2_noxjk-z2Eq2_noxjk))+1);
                            P_raEq2_noxjk=P_raEq2_noxjk+(pr_action2Eq2_noxjk/odors);
                        end    
                        
                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1Eq2_wHy= 1/(1+exp(C*(z2Eq2_wHy-z1Eq2_wHy)));
                            
                        
                            p_raEq2_wHy=p_raEq2_wHy+(pr_action1Eq2_wHy/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2Eq2_wHy=  1/(exp(C*(z1Eq2_wHy-z2Eq2_wHy))+1);
                            p_raEq2_wHy=p_raEq2_wHy+(pr_action2Eq2_wHy/odors);
                        end 
                        
                         
                        
                    end
end



end