%Calculation of Stationary Distribution

function [st,FOUND,max_diff,error] = Stationary(E,iteration_num,st_old)

    display('Inside Stationary');
    TRUE=1;
    FALSE=0;
    [V lambda]=eig(E');
    %display(E);
    %display(lambda);
    %display(V);
     
    %search for 1 in lambda
    row_prime=-1;
    for row=1:3
        if( int32(lambda(row,row))==1) %as lambda is a diagonal matri
        %if((10000*lambda(row,row)) ==10000)
           row_prime=row;
           break;
        end      
    end
           
    st=V(:,row_prime);
    display(st);
    sum_eigen=0;
    %check if at is a mixture of positive and negative numbers
    flag=0;
    error=0;
    for row=1:3
       if(flag==0 && st(row)> 0 )
          flag=1;
       elseif(flag==0 && st(row) < 0 )
          flag=-1;
       elseif(flag==1 && st(row) < 0 )
          error=1;
       elseif(flag==-1  && st(row) > 0 )
          error=1;
       end
       
       sum_eigen=sum_eigen+st(row);
   end%end for row

    if(error==1) 
       disp('iteration'); disp(iteration_num); disp(' (error');
       disp(E);
    end
           
   for row=1:3
       st(row)=st(row)/sum_eigen;
   end

      
  FOUND=TRUE;
   if(iteration_num==1) % exception for iteration 1
      FOUND=FALSE;
   else
      for row=1:3
          if(abs(st(row)- st_old(row))>0.00001) 
               FOUND=FALSE;
          end
      end
   end

   if(FOUND==TRUE) 
      text=sprintf('  Ending in iteration number %d',iteration_num);
 	  disp(text);
   end
               
   max_val=-1;
   min_val=99;
   for row=1:3
       if( max_val<=st(row))
           max_val=st(row);
       end
       if ( min_val>=st(row))
            min_val=st(row);
        end
   end
   max_diff=max_val-min_val;
   text=sprintf(' %d) max_diff=%d',iteration_num,max_diff);
   disp(text);
end %end of function
