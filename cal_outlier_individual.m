function [num_val] = cal_outlier_individual(x,x_5,j,row,NO_AREA,deg,num)
%num = index of the neighbour
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %         find the top f and bottom f outliers
       %         max_sum and min_sum are the top and bottom outliers
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       sum_violation=zeros(1,NO_AREA-1); %NO_AREA-1 because one does not use one's own row
       num_val=zeros(deg,1);
       
       for i=1:deg
           pos=1;
           for col=1:NO_AREA
               if(col~=row)
                    if(col~=5)
                        sum_temp=abs(x((row-1)*deg+i,j)-x((col-1)*deg+i,j));
                    else
                        sum_temp=abs(x((row-1)*deg+i,j)-reshape(x_5(row,i,j),1,1));
                    end
               
                    sum_violation(1,pos)=sum_temp;
                    pos=pos+1;
               end %end if(col~=row)
           end
           [~, max_index ]=max(sum_violation);
           [~, min_index ]=min(sum_violation);
           if(num==max_index)
                num_val(i,1)=x((row-1)*deg+i,j);
           elseif (num==min_index)
                num_val(i,1)=x((row-1)*deg+i,j);
           elseif(num~=5)
                num_val(i,1)=x((num-1)*deg+i,j);
           else %num==5
                num_val(i,1)=reshape(x_5(row,i,j),1,1);
           end
       end %end for i=1 to deg
        
