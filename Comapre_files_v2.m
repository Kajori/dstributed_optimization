clear all; clc; close all;

TOLERANCE=0;

load('fault_detection_3per.mat');
load('fault_no_detection_3per.mat');


var_v9=x_fault_no_detection;
var_v10=x_fault_detection;
temp=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compare Version area 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dummy=abs(var_v9(1,1)-var_v10(1,1));
max_error_percentage= (dummy/abs(var_v10(1,1)));  
display(max_error_percentage);
N=10000;

col=N;
for row=1:size(var_v9,1)
      
        dummy=abs(var_v9(row,col)-var_v10(row,col));
        if(dummy~=0)
            text=sprintf('\n dummy=%d',dummy);
                 disp(text);
          error_percentage= (dummy/abs(var_v10(row,col)));
          if(max_error_percentage<=error_percentage) %;%temp_8=abs(temp_9-Z)<=TOLERANCE;%temp_8=isequal(x(1:1*deg,j+1),Z(:));
                 temp=0;
                  
                 error_col=row-(floor(row/deg)*deg);
                 text=sprintf('\n ROW=%d \n Error row=%d col= %d \n diff=%d',row,floor(row/deg+1),error_col,error_percentage);
                 disp(text);
                 text=sprintf(' Original values \n%d \n %d ',var_v9(row,col),var_v10(row,col));
                 disp(text);
                 
                 max_error_percentage=error_percentage;
          end
        end
   
end

display(temp);
display(max_error_percentage);

  plot(0:size(var_v10,2)-1,var_v10(39,1:end),'linewidth',2);
  hold on
  plot(0:size(var_v9,2)-1,var_v9(39,1:end),'linewidth',2);
  
 
