function [p_raH]=HomogenousModel_KernelTesting(C,WopHom,PNtrials, thetaH_Ftheta,classAction1,numTrials,numtrainingSamples,YHomog)



numm=numTrials-numtrainingSamples;

n=2000;

p_raH=0;



odors= (size(YHomog,2))*numm;



for trials=numtrainingSamples+1:numTrials

                    for odour=1:(size(YHomog,2))

                        

                        z1H=WopHom(:,1)'*YHomog(:,odour,trials);
                        z2H=WopHom(:,2)'*YHomog(:,odour,trials);
                        
                        %% prob. of success in the homogenous model 
                        if ( ~ isempty(find(classAction1==odour)) )
                            
                            pr_action1H= 1/(1+exp(C*(z2H-z1H)));
                            
                            %if(pr_action1H>rand(1))
                            p_raH=p_raH+(pr_action1H/odors);
                            
                            
                        end

                        if (  isempty((find(classAction1==odour))) )
                            
                            pr_action2H=  1/(exp(C*(z1H-z2H))+1);
                            %if(pr_action2H>rand(1))
                                p_raH=p_raH+(pr_action2H/odors);
                            
                            
                        end
                        
                       
                        
                        
                        
                        
                        
                        
                        
                    end
end



end