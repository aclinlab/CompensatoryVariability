function [p_ra,p_raEq,p_raEq2,p_raThetaHomeo, p_raInhPlast,p_raEq2_noxjk]=Robustness_testing(C,Wop,WopEq,WopEq2,WopThetaHoemo, WopInhPlast,WopEq2_noxjk,...
    classAction1,numTrials,numtrainingSamples,Y,YEqualized,Y_comp2, Y_theta_homeo, Y_inhPlast,Y_comp2_noxjk)

              

numm=numTrials-numtrainingSamples;
    

p_ra=0;
p_raEq=0;
p_raEq2=0;
p_raInhPlast=0;
p_raThetaHomeo=0;
p_raEq2_noxjk=0;

% p_raEq2_wHY=0;
% p_raEq2_noxjk=0;
% 
% p_raKenn_CL_HY=0;
% p_raKenn_CL_noHY=0;
% 
% p_raInhPlast_wHY=0;


i=1; 

odors= (size(Y,2))*numm;

for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(Y,2))

                        z1=Wop(:,1)'*Y(:,odour,trials);
                        z2=Wop(:,2)'*Y(:,odour,trials);

                        z1Eq=WopEq(:,1)'*YEqualized(:,odour,trials);
                        z2Eq=WopEq(:,2)'*YEqualized(:,odour,trials);

                        %%
                        z1Eq2=WopEq2(:,1)'*Y_comp2(:,odour,trials);
                        z2Eq2=WopEq2(:,2)'*Y_comp2(:,odour,trials);
                        
                        z1Eq2_noxjk=WopEq2_noxjk(:,1)'*Y_comp2_noxjk(:,odour,trials);
                        z2Eq2_noxjk=WopEq2_noxjk(:,2)'*Y_comp2_noxjk(:,odour,trials);
                        
%                         z1Eq2_wHY=WopEq2_HY(:,1)'*Y_comp2_wHY(:,odour,trials);
%                         z2Eq2_wHY=WopEq2_HY(:,2)'*Y_comp2_wHY(:,odour,trials);
%                         
%                         z1Eq2_noxjk=WopEq2_noxjk(:,1)'*Y_comp2_noxjk(:,odour,trials);
%                         z2Eq2_noxjk=WopEq2_noxjk(:,2)'*Y_comp2_noxjk(:,odour,trials);

                        %%                       
                       
                        
                        z1ThetaHome=WopThetaHoemo(:,1)'*Y_theta_homeo(:,odour,trials);
                        z2ThetaHome=WopThetaHoemo(:,2)'*Y_theta_homeo(:,odour,trials);
                        
%                         z1Kenn_CL_HY=WopKenn_CL_HY(:,1)'*Y_Kenn_CL_HY(:,odour,trials);
%                         z2Kenn_CL_HY=WopKenn_CL_HY(:,2)'*Y_Kenn_CL_HY(:,odour,trials);
%                         
%                         z1Kenn_CL_noHY=WopKenn_CL_noHY(:,1)'*Y_Kenn_CL_noHY(:,odour,trials);
%                         z2Kenn_CL_noHY=WopKenn_CL_noHY(:,2)'*Y_Kenn_CL_noHY(:,odour,trials);
                        %%
                        z1InhPlast=WopInhPlast(:,1)'*Y_inhPlast(:,odour,trials);
                        z2InhPlast=WopInhPlast(:,2)'*Y_inhPlast(:,odour,trials);
%                         
%                         z1InhPlast_wHY=WopInhPlast_HY(:,1)'*Y_inhibPlast_HY(:,odour,trials);
%                         z2InhPlast_wHY=WopInhPlast_HY(:,2)'*Y_inhibPlast_HY(:,odour,trials);
%                         
                        %%
                        
                        if (~ isempty(find(classAction1==odour)) )
                            
                            
                            pr_action1=  1/(1+exp(C*(z2-z1)));
                            
                            p_ra=p_ra+(pr_action1/odors);
                            %% store the average difference in activity for every right action
                            
                            avgdA(i)= abs(z2-z1);
                            i=i+1;
                            
                        end

                        if ( isempty((find(classAction1==odour))) )
                            
                            
                            pr_action2=  1/(exp(C*(z1-z2))+1);
                            
                            
                          
                            p_ra=p_ra+(pr_action2/odors);
                            
                            avgdA(i)= abs(z2-z1);
                            i=i+1;
                        
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
                         
                        
                        %    
%                         if ( ~ isempty(find(classAction1==odour)) )
%                             
%                             pr_action1Eq2_wHY= 1/(1+exp(C*(z2Eq2_wHY-z1Eq2_wHY)));
%                             
%                         
%                             p_raEq2_wHY=p_raEq2_wHY+(pr_action1Eq2_wHY/odors);
%                             
%                             
%                         end
% 
%                         if (  isempty((find(classAction1==odour))) )
%                             
%                             pr_action2Eq2_wHY=  1/(exp(C*(z1Eq2_wHY-z2Eq2_wHY))+1);
%                             p_raEq2_wHY=p_raEq2_wHY+(pr_action2Eq2_wHY/odors);
%                                
%                                 
%                         end
%                         
%                         %
%                         %    
%                         if ( ~ isempty(find(classAction1==odour)) )
%                             
%                             pr_action1Eq2_noxjk= 1/(1+exp(C*(z2Eq2_noxjk-z1Eq2_noxjk)));
%                             
%                         
%                             p_raEq2_noxjk=p_raEq2_noxjk+(pr_action1Eq2_noxjk/odors);
%                             
%                             
%                         end
% 
%                         if (  isempty((find(classAction1==odour))) )
%                             
%                             pr_action2Eq2_noxjk=  1/(exp(C*(z1Eq2_noxjk-z2Eq2_noxjk))+1);
%                             p_raEq2_noxjk=p_raEq2_noxjk+(pr_action2Eq2_noxjk/odors);
%                                
%                                 
%                         end
                        
                        %
                        
                        
                        
                         if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1thetahomeo= 1/(1+exp(C*(z2ThetaHome-z1ThetaHome)));
                            
                            p_raThetaHomeo=p_raThetaHomeo+(pr_action1thetahomeo/odors);
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2thetahomeo=  1/(exp(C*(z1ThetaHome-z2ThetaHome))+1);
                            p_raThetaHomeo=p_raThetaHomeo+(pr_action2thetahomeo/odors);      
                        end
                        
                        
                         if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1Eq2_noxjk= 1/(1+exp(C*(z2Eq2_noxjk-z1Eq2_noxjk)));
                            
                        
                            p_raEq2_noxjk=p_raEq2_noxjk+(pr_action1Eq2_noxjk/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2Eq2_noxjk=  1/(exp(C*(z1Eq2_noxjk-z2Eq2_noxjk))+1);
                            p_raEq2_noxjk=p_raEq2_noxjk+(pr_action2Eq2_noxjk/odors);
                        end 
                        
                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1InhPlast= 1/(1+exp(C*(z2InhPlast-z1InhPlast)));
                            
                        
                            p_raInhPlast=p_raInhPlast+(pr_action1InhPlast/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2InhPlast=  1/(exp(C*(z1InhPlast-z2InhPlast))+1);
                            p_raInhPlast=p_raInhPlast+(pr_action2InhPlast/odors);
                                
                               
                                
                        end
                        
%                         if ( ~ isempty(find(classAction1==odour)) )
%                             
%                             pr_action1InhPlast_wHY= 1/(1+exp(C*(z2InhPlast_wHY-z1InhPlast_wHY)));
%                             
%                         
%                             p_raInhPlast_wHY=p_raInhPlast_wHY+(pr_action1InhPlast_wHY/odors);
%                             
%                             
%                         end
% 
%                         if (  isempty((find(classAction1==odour))) )
%                             
%                             pr_action2InhPlast_wHY=  1/(exp(C*(z1InhPlast_wHY-z2InhPlast_wHY))+1);
%                             p_raInhPlast_wHY=p_raInhPlast_wHY+(pr_action2InhPlast_wHY/odors);
%                                 
%                             
%                                 
%                                 
%                         end
                        
                    end
end



end