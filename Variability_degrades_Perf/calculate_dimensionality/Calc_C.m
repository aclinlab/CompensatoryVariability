function C= Calc_C (h)
 M=size(h,1);
 h_avg= mean(h,2);
 C=zeros(M,M);
 h_avgMat=repmat(h_avg,1,size(h,2));
 
%for i=1:M
 %  for j=1:M
      
C= (h-h_avgMat)*(h-h_avgMat)'; 
 
%      hi= h(i,:,:); % 1xodsxnum trials.  
%      hi_minus_hiavg= hi-h_avg(i);
%      
%      hj= h(j,:,:); % 1xodsxnum trials. 
%      hj_minus_hjavg= hj-h_avg(j);
%      
%      temp=hi_minus_hiavg .* hj_minus_hjavg;
%      C(i,j)=mean(temp(:));
       
  % end
    
%end




end