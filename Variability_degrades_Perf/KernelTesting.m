function [p_ra,p_raH,pr_aVarW,pr_aVarN]=KernelTesting(C,Wop,WopHom, WopVarW_fixedN,WopVarN_fixedW , classAction1,numTrials,numtrainingSamples,Y,YHomog, Y_varw_fixedN, Y_varN_fixedW)


numm=numTrials-numtrainingSamples;
n=2000;
p_ra=0;
p_raH=0;
pr_aVarW=0;
pr_aVarN=0;
odors= (size(Y,2))*numm;

for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(Y,2))

                        z1=Wop(:,1)'*Y(:,odour,trials);
                        z2=Wop(:,2)'*Y(:,odour,trials);

                        z1H=WopHom(:,1)'*YHomog(:,odour,trials);
                        z2H=WopHom(:,2)'*YHomog(:,odour,trials);
                        
                        z1VarW= WopVarW_fixedN(:,1)'*Y_varw_fixedN(:,odour,trials);
                        z2VarW= WopVarW_fixedN(:,2)'*Y_varw_fixedN(:,odour,trials);
                        
                        z1VarN=WopVarN_fixedW(:,1)'*Y_varN_fixedW(:,odour,trials);
                        z2VarN=WopVarN_fixedW(:,2)'*Y_varN_fixedW(:,odour,trials);
                        

                        if (~ isempty(find(classAction1==odour)) )
                                                        
                            pr_action1=  1/(1+exp(C*(z2-z1)));                            
                            p_ra=p_ra+(pr_action1/odors);
                                                      
                        end

                        if ( isempty((find(classAction1==odour))) )
                            
                            pr_action2=  1/(exp(C*(z1-z2))+1);                                                                                  
                            p_ra=p_ra+(pr_action2/odors);
                            
                        end
                        
                        %% prob. of success in the homogenous model 
                        if ( ~ isempty(find(classAction1==odour)) )
                            pr_action1H= 1/(1+exp(C*(z2H-z1H)));                            
                            p_raH=p_raH+(pr_action1H/odors);
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2H=  1/(exp(C*(z1H-z2H))+1);
                            p_raH=p_raH+(pr_action2H/odors);
                           
                        end
                        
                        %% prob. of success in the var w fixed N model
                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1VarW= 1/(1+exp(C*(z2VarW-z1VarW)));
                            
                        
                            pr_aVarW=pr_aVarW+(pr_action1VarW/odors);
  
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2VarW=  1/(exp(C*(z1VarW-z2VarW))+1);
                            
                            pr_aVarW=pr_aVarW+(pr_action2VarW/odors);
     
                        end
                        
                         %% prob. of success in the var N fixed W model
                        
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1VarN= 1/(1+exp(C*(z2VarN-z1VarN)));
                            
                        
                            pr_aVarN=pr_aVarN+(pr_action1VarN/odors);
  
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2VarN=  1/(exp(C*(z1VarN-z2VarN))+1);
                            
                            pr_aVarN=pr_aVarN+(pr_action2VarN/odors);
     
                        end
                       
                    end
end

end