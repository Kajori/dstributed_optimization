function [mal_node_1,mal_node_2,Pre_succ] = cal_outlier_median_v1(x,j,row,NO_AREA,deg)

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %         calculate the median
       %         remove the ouliers from the median
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       med=zeros(1,deg);
       
       temp=zeros(NO_AREA,deg);
       for col=1:deg
           for row_temp=1:NO_AREA
                temp(row_temp,col)=x((row_temp-1)*deg+col);
           end
       end
       M = median(temp);
       if(j==1) disp('size(M)');  disp(size(M)); disp('size(temp)'); disp(size(temp)); end
       
       sum_violation=zeros(1,NO_AREA);
       for col=1:NO_AREA
            sum_temp=0.5*norm(M'-x((col-1)*deg+1:col*deg,j));
            sum_violation(1,col)=sum_temp;
       end
       [max_val, max_index ]=max(sum_violation);%added Rakesh's suggestion
       [min_val, min_index ]=min(sum_violation); %added Rakesh's suggestion
       %check for the presence of another max_val
       
       count=0;
       for col=1:NO_AREA
           if (max_val==sum_violation(1,col)) count=count+1; end
       end
       if (count>1) max_index=-1; end
       
       count=0;
       for col=1:NO_AREA
           if (min_val==sum_violation(1,col)) count=count+1; end
       end
       if (count>1) min_index=-1; end
           
           
       mal_node_1=-1;
       mal_node_2=-1;

       if (row==1 && max_index==5) mal_node_1=5;  end
       if (row==1 && min_index==5) mal_node_1=5;  end
       if (row==1 && max_index==2) mal_node_2=2;  end
       if (row==1 && min_index==2) mal_node_2=2;  end
       
       if (row==2 && max_index==1) mal_node_1=1;  end
       if (row==2 && min_index==1) mal_node_1=1;  end
       if (row==2 && max_index==3) mal_node_2=3;  end
       if (row==2 && min_index==3) mal_node_2=3;  end
       
       if (row==3 && max_index==2) mal_node_1=2;  end
       if (row==3 && min_index==2) mal_node_1=2;  end
       if (row==3 && max_index==4) mal_node_2=4;  end
       if (row==3 && min_index==4) mal_node_2=4;  end
           
           
       if (row==4 && max_index==3) mal_node_1=3;  end
       if (row==4 && min_index==3) mal_node_1=3;  end
       if (row==4 && max_index==5) mal_node_2=5;  end
       if (row==4 && min_index==5) mal_node_2=5;  end
       
       if (row==5 && max_index==4) mal_node_1=4;  end
       if (row==5 && min_index==4) mal_node_1=4;  end
       if (row==5 && max_index==1) mal_node_2=1;  end
       if (row==5 && min_index==1) mal_node_2=1;  end
         
       Pre_succ=2;
       if(row==1 && mal_node_1==5)     Pre_succ=Pre_succ-1; end
       if(row==1 && mal_node_2==2)     Pre_succ=Pre_succ-1; end
       
       if(row==2 && mal_node_1==1)     Pre_succ=Pre_succ-1; end
       if(row==2 && mal_node_2==3)     Pre_succ=Pre_succ-1; end
       
       if(row==3 && mal_node_1==2)     Pre_succ=Pre_succ-1; end
       if(row==3 && mal_node_2==4)     Pre_succ=Pre_succ-1; end
       
       if(row==4 && mal_node_1==3)     Pre_succ=Pre_succ-1; end
       if(row==4 && mal_node_2==5)     Pre_succ=Pre_succ-1; end
       
       if(row==5 && mal_node_1==4)     Pre_succ=Pre_succ-1; end
       if(row==5 && mal_node_2==1)     Pre_succ=Pre_succ-1; end
