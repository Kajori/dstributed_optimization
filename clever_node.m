% there if the neighbours are malacious its ets the node 1 and node 2
% respectively fot the left and right neighbours
%xsent_1 is the value sent by 5 to neighbour 1
function [xsent_1,xsent_4] = clever_node(my_row,x,j,deg)

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %         find the top f and bottom f outliers
       %         max_sum and min_sum are the top and bottom outliers
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       sum_violation=zeros(1,NO_AREA-1);
       row=1;
       for col=1:NO_AREA-1 %we wil not condsider 5
            sum_temp=0.5*norm(x((row-1)*deg+1:row*deg,j)-x((col-1)*deg+1:col*deg,j));
            sum_violation(1,col)=sum_temp;
       end
       [max_val,~ ]=max(sum_violation);
       
       temp=x((row-1)*deg+1:row*deg,j+1);
       
    
