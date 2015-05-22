% there if the neighbours are malacious its ets the node 1 and node 2
% respectively fot the left and right neighbours
%xsent_1 is the value sent by 5 to neighbour 1
% it is written for dist_markov_Chain_imp3
function [x_5] = clever_adversary(my_row,x_5,NO_AREA,x,j,deg)
     disp(sprintf(' -------------- CA %d --------------',j));
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %         find the top f and bottom f outliers
       %         max_sum and min_sum are the top and bottom outliers
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       for row=1:NO_AREA % we need to cal x5 for each node
           if(row~=my_row)
               for i=1:deg % for each element
                    sum_violation=zeros(1,NO_AREA-2);
                    pos=1;
                    for col=1:NO_AREA-1 %compare two elements %we do not do it for the row 5
                        if(col~=row)
                            sum_temp=abs(x((row-1)*deg+i,j)-x((col-1)*deg+i,j));
                            sum_violation(1,pos)=sum_temp;
                            pos=pos+1;
                        end %end if(col~=row)
                    end %for col
                    [max_val, max_index ]=max(sum_violation);
                    [min_val, min_index ]=min(sum_violation);
                    disp('min_val');
                    disp(min_val);
                    disp('sum_violation');
                    disp(sum_violation);
                    temp=((max_val-min_val)*randn(1))+min_val;
                    disp('temp');
                    disp(temp);
                    x_5(row,i,j)=temp;
               end%for i deg
           end %if(row~=my_row)
       end %for row
        

       
    
