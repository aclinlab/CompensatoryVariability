function [p_ra,p_raEq,p_raEq2,p_raKenn,p_raInhPlast,p_raThetaHomeo]=testingModels_accuracies(C,Wop,WopEq,WopEq2, WopKenn , WopInhPlast,WopThetaHoemo,...
    classAction1,numTrials,numtrainingSamples,Y,YEqualized,Y_comp2,Y_Kenn,Y_inhPlast, Y_theta_homeo)

              
numm=numTrials-numtrainingSamples;
p_ra=0;
p_raEq=0;
p_raEq2=0;
p_raInhPlast=0;
p_raThetaHomeo=0;
p_raKenn=0;

i=1; 

odors= (size(Y,2))*numm;

for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(Y,2))

                        z1=Wop(:,1)'*Y(:,odour,trials);
                        z2=Wop(:,2)'*Y(:,odour,trials);

                        z1Eq=WopEq(:,1)'*YEqualized(:,odour,trials);
                        z2Eq=WopEq(:,2)'*YEqualized(:,odour,trials);
                     
                        z1Eq2=WopEq2(:,1)'*Y_comp2(:,odour,trials);
                        z2Eq2=WopEq2(:,2)'*Y_comp2(:,odour,trials);
                        
%                         z1Eq2_noxjk=WopEq2_noxjk(:,1)'*Y_comp2_noxjk(:,odour,trials);
%                         z2Eq2_noxjk=WopEq2_noxjk(:,2)'*Y_comp2_noxjk(:,odour,trials);                        

                        z1ThetaHome=WopThetaHoemo(:,1)'*Y_theta_homeo(:,odour,trials);
                        z2ThetaHome=WopThetaHoemo(:,2)'*Y_theta_homeo(:,odour,trials);
                        
                        z1Kenn=WopKenn(:,1)'*Y_Kenn(:,odour,trials);
                        z2Kenn=WopKenn(:,2)'*Y_Kenn(:,odour,trials);                        
                      
                        z1InhPlast=WopInhPlast(:,1)'*Y_inhPlast(:,odour,trials);
                        z2InhPlast=WopInhPlast(:,2)'*Y_inhPlast(:,odour,trials);
            
                        %%
                        
                        if (~ isempty(find(classAction1==odour)) )
                            
                            
                            pr_action1=  1/(1+exp(C*(z2-z1)));
                            
                            p_ra=p_ra+(pr_action1/odors);

                        end

                        if ( isempty((find(classAction1==odour))) )
                            
                            
                            pr_action2=  1/(exp(C*(z1-z2))+1);
                            
            
                            p_ra=p_ra+(pr_action2/odors);

                        end


                        %% prob. of success in equalized model
                        if ( ~ isempty(find(classAction1==odour)) )
                           
                            pr_action1Eq=  1/(1+exp(C*(z2Eq-z1Eq)));
                            
                            %if(pr_action1Eq>rand(1))
                            p_raEq=p_raEq+(pr_action1Eq/odors);
                                                           
                        end

                        if ( isempty((find(classAction1==odour))) )
                            
                            pr_action2Eq=  1/(exp(C*(z1Eq-z2Eq))+1);
                            
                           p_raEq=p_raEq+(pr_action2Eq/odors);
                                                            
                        end
                        
                       
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1Eq2= 1/(1+exp(C*(z2Eq2-z1Eq2)));
                            
                        
                            p_raEq2=p_raEq2+(pr_action1Eq2/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2Eq2=  1/(exp(C*(z1Eq2-z2Eq2))+1);
                            p_raEq2=p_raEq2+(pr_action2Eq2/odors);
                        end     
                         

                         if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1Kenn= 1/(1+exp(C*(z2Kenn-z1Kenn)));
                            
                            p_raKenn=p_raKenn+(pr_action1Kenn/odors);
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2Kenn=  1/(exp(C*(z1Kenn-z2Kenn))+1);
                            p_raKenn=p_raKenn+(pr_action2Kenn/odors);      
                        end
                        
                        
                        
                         if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1thetahomeo= 1/(1+exp(C*(z2ThetaHome-z1ThetaHome)));
                            
                            p_raThetaHomeo=p_raThetaHomeo+(pr_action1thetahomeo/odors);
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2thetahomeo=  1/(exp(C*(z1ThetaHome-z2ThetaHome))+1);
                            p_raThetaHomeo=p_raThetaHomeo+(pr_action2thetahomeo/odors);      
                        end
                        

                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1InhPlast= 1/(1+exp(C*(z2InhPlast-z1InhPlast)));
                            
                        
                            p_raInhPlast=p_raInhPlast+(pr_action1InhPlast/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2InhPlast=  1/(exp(C*(z1InhPlast-z2InhPlast))+1);
                            p_raInhPlast=p_raInhPlast+(pr_action2InhPlast/odors);
                                
                               
                                
                        end
                        
                      
                        
                    end
end



end