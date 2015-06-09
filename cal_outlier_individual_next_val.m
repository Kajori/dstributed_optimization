%here the malacious node sends the same value to all the other nodes
%Remove min-max of each eigen value and replace with next highest/lowest
function [num_val] = cal_outlier_individual_next_val(x,j,row,NO_AREA,deg,num,PRECISION)
%num = index of the neighbour
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %         find the top f and bottom f outliers
       %         max_sum and min_sum are the top and bottom outliers
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       sum_violation=zeros(1,NO_AREA); %NO_AREA because one does not use one's own row
       num_val=zeros(deg,1);
       
       for i=1:deg
           for col=1:NO_AREA
                    sum_temp=x((col-1)*deg+i,j);
                    sum_violation(1,col)=sum_temp;
           end %end of col=1:NO_AREA
           
          
           [max_val, max_index ]=max(sum_violation);
           [min_val, min_index ]=min(sum_violation);
           

           if(num==max_index)
               count=0;
               for col=1:NO_AREA
                   if(max_val==x((col-1)*deg+i,j))
                       count=count+1;
                   end
               end
               if(count>1)
                   num_val(i)=x((num-1)*deg+i,j);
               else
                   sum_violation(num)=-9999;
                   [next_max_val,~]=max(sum_violation);
                   num_val(i)=next_max_val;
               end
           elseif(num==min_index)  
               count=0;
               for col=1:NO_AREA %Replace with next highest/lowest
                   if(min_val==x((col-1)*deg+i,j))
                       count=count+1;
                   end
               end
               if(count>1)
                   num_val(i)=x((num-1)*deg+i,j);
               else %Replace with next highest/lowest
                   sum_violation(num)=9999;
                   [next_min_val,~]=min(sum_violation);
                   num_val(i)=next_min_val;
               end
           else
                num_val(i)=x((num-1)*deg+i,j);
           end
           
           if (j>10000)   disp(sum_violation); disp(sprintf(' j= %d i=%d row=%d num=%d max_index=%d min_index=%d',j,i,row,num,max_index,min_index)); end
     end %end for i=1 to deg
        
