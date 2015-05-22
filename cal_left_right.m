function [left,right] = cal_left_right(row)
      %find the left and right indexes
      %it fllows a ring fashion
      
       if(row>=2 && row<=4)
           left=row-1;
           right=row+1;
       elseif(row==1)
           left=5;
           right=2;
       elseif (row==5)
           left=4;
           right=1;
       end

      
