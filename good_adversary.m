% there if the neighbours are malacious its ets the node 1 and node 2
% respectively fot the left and right neighbours
% here the adversary does not manipulate nay values
% it is written for dist_markov_Chain_imp3
function [x_5] = good_adversary(my_row,x_5,NO_AREA,x,j,deg)
       for row=1:NO_AREA % we to to cal x5 for each node
           if(row~=my_row)
               for i=1:deg % for each element
                   x_5(row,i,j)=x((my_row-1)*deg+i,j);
               end%for i deg
           end %if(row~=my_row)
       end %for row
        

       
    
