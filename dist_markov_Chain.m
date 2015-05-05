%3 phase fault  -no
%Gaussian Noise - YES
% malacious node - if present node 5
% there will be a central guy.. trim the top f and bottom f values 
%here f =1
MAX_DEV = 1e-6;

clear all; clc; close all;
%load n16 % PST model  
%load IEEE68_dyn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Defining the constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NO_AREA=5;
PMU_PER_AREA=3;
SIZE_OF_voltage_data=1000;
NO_OF_EIGEN_VALUES=20;

PRECISION=10000000000;
TRUE=1;
FALSE=0;
st_prev=[0,0,0];
area_line_no=[53,60,61;
              30,48,63;
              62,66,67;
              64,65,68;
              54,56,57];
no_of_successor  =[2,1,1,1,0];
no_of_predecessor=[0,1,1,1,2];

successor=[2,5,-1;
      3,-1,-1;
      4,-1,-1;
      5,-1,-1];
predecessor=[-1,-1;
              1,-1;
              2,-1;
              3,-1;
              1,4];
          
Y=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);
Y_org=zeros(NO_AREA,PMU_PER_AREA,SIZE_OF_voltage_data);
deg=2*NO_OF_EIGEN_VALUES; 
m=60; % height of the Hankel matrix
H_pmu=zeros(NO_AREA,PMU_PER_AREA,m,deg);
C_pmu=zeros(NO_AREA,PMU_PER_AREA,m);

H_pdc=zeros(NO_AREA,PMU_PER_AREA*m,deg);
C_pdc=zeros(NO_AREA,PMU_PER_AREA*m);
display('Gaussian Noise no 3 phase fault');







rho=10^-9;
N=10000;
disp('-----------------------------------------------');
tic;

TT=PMU_PER_AREA*m;
mal_node=-1;
FOUND=0;
SCALE = 0.001; %0.01 % noise
E1=zeros(3,3);E2=zeros(3,3);E3=zeros(3,3);E4=zeros(3,3);E5=zeros(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Loading data from the input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
for row=1:NO_AREA
     for col=1:PMU_PER_AREA
           file_name=sprintf('no_3_phase_fault_area_%d_line_%d.txt',row,area_line_no(row,col));
           Y_org(row,col,:)=importdata(file_name,'\n');
      end %end of for col
end %end of for row
disp(' Loaded data from the input files ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Write statistics to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('stat.txt', 'a+');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %         adding Gausian Noise  to the data
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row=1:NO_AREA 
    for col=1:PMU_PER_AREA
             noise = randn(1); % noise with mean=0 and std=1;
             Y(row,col,:)= ( 1+  noise* SCALE )* Y_org(row,col,:);
    end %end of for col
end %end of for row        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Hankel Matrix Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row=1:NO_AREA
    for col=1:PMU_PER_AREA
               temp_1=reshape(Y(row,col,:),SIZE_OF_voltage_data,1);
               temp_2=fliplr(Hankel(temp_1,m,deg));
               H_pmu(row,col,:,:)=temp_2;
               temp_3=reshape(Y(row,col,1:m),1,m);
               C_pmu(row,col,:)=temp_3.';      
    end %end of for col
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Combining Data from PMU to PDC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for col=1:PMU_PER_AREA
          H_pdc(row,(col-1)*m+1: col*m,:)= reshape(H_pmu(row,col,:,:),m,deg);
          C_pdc(row,(col-1)*m+1: col*m) = reshape(C_pmu(row,col,:),m,1);
     end %end of for col
 end %end of for row

w=zeros(NO_AREA*deg,1);
x(:,1)=ones(NO_AREA*deg,1); %Everyone has equal initial x

for j=1:N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Check if we have reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( convergence_flag==TRUE)
         text=sprintf('Convergence Reached at %d ',j);
         display(text);
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         if we have not reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text=sprintf(' -------------- %d --------------',j);
    disp( text);
    
    H_r=reshape(H_pdc(row,1:TT,:),TT,deg);
    C_r=reshape(C_pdc(row,1:TT,:),TT,deg);
    
    x((row-1)*deg+1:row*deg,j+1)=(H_r'*H_r+rho*eye(size(H_r,2)))\(H_r'*C_r-w((row-1)*deg+1:row*deg,j)+rho*z(:,j));
    
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
   
     Pre_succ=2; %no_of_successor(row)+no_of_predecessor(row); % Pre_succ=n(Pi)+n(Si)
       
     dummy=reshape(H_pdc(row,1:TT,:),TT,deg);
     dummy_0=dummy.';
     dummy_01=dummy_0*dummy;
        
     if(row==1 && mal_node==5)
         Pre_succ=Pre_succ-1;
     elseif(row==4 && mal_node==5)
          Pre_succ=Pre_succ-1;
     end
     
     dummy_1=dummy_01+rho*Pre_succ*eye(size(dummy,2));
     dummy_3=reshape(C_pdc(row,1:TT,:),TT,1);
     dummy_4=dummy'*dummy_3;
     mal_node=-1;  %for the global state calculation %REMOVE FOR CHECKING
         
     if(row==1)  % we have enabled chaking oonly for node 2 i.e. only node 5 can be malacious
       if  (mal_node~=5)
            dummy_5=dummy_4+(w12(:,j) + w15(:,j))+rho*(x(deg+1:2*deg,j)+x(4*deg+1:5*deg,j));
        else
             dummy_5=dummy_4+(w12(:,j))+rho*(x(deg+1:2*deg,j));
        end
        x((row-1)*deg+1:row*deg,j+1) =dummy_1\dummy_5;
     
     elseif(row==2)
        dummy_5=dummy_4+(w23(:,j) - w12(:,j))+rho*(x(1:deg,j+1)+x(2*deg+1:3*deg,j));
        x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
        w12(:,j+1)=w12(:,j)-rho*(x(1:deg,j+1)-x(deg+1:2*deg,j+1));   % w12 update as 2 is the successor of 1
     
     elseif(row==3)
         dummy_5=dummy_4+(w34(:,j) - w23(:,j))+rho*(x(deg+1:2*deg,j+1)+x(3*deg+1:4*deg,j));
         x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
         w23(:,j+1)=w23(:,j)-rho*(x(deg+1:2*deg,j+1)-x(2*deg+1:3*deg,j+1));   % w23 update as 2 is the predecessor of 3
     
     elseif(row==4)
          if (mal_node~=5)
              dummy_5=dummy_4+(w45(:,j) - w34(:,j))+rho*(x(2*deg+1:3*deg,j+1)+x(4*deg+1:5*deg,j));
          else
              dummy_5=dummy_4+( - w34(:,j))+rho*(x(2*deg+1:3*deg,j));
          end
          x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
          w34(:,j+1)=w34(:,j)-rho*(x(2*deg+1:3*deg,j+1)-x(3*deg+1:4*deg,j+1));        % w34 update as 3 is the predecessor of 4     
            
    elseif ( row==5  && row~=mal_node)
           dummy_5=dummy_4+(-w15(:,j) - w45(:,j))+rho*(x(1:deg,j+1)+x(3*deg+1:4*deg,j+1));
           x((row-1)*deg+1:row*deg,j+1)=1.05*dummy_1\dummy_5;
           %x((row-1)*deg+1:row*deg,j+1)=1.05*  x((row-1)*deg+1:row*deg,j+1);
           w15(:,j+1)=w15(:,j)-rho*(x(1:deg,j+1)-x(4*deg+1:5*deg,j+1));        
           w45(:,j+1)=w45(:,j)-rho*(x(3*deg+1:4*deg,j+1)-x(4*deg+1:5*deg,j+1));      
        end %end of if
   end %end of for row
   
    convergence_flag=TRUE;
    for row=1:NO_AREA*deg
        if ((abs(x(row,j+1)-x(row,j)))>0.000001)
          convergence_flag=FALSE;
        end
    end
    if( convergence_flag==TRUE)
         display('2. convergence_flag==TRUE');
    end 
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %         after receiving from all successor
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w12(:,j+1)=w12(:,j)-rho*(x(1:deg,j+1)-x(deg+1:2*deg,j+1));
    if(mal_node~=5)
          w15(:,j+1)=w15(:,j)-rho*(x(1:deg,j+1)-x(4*deg+1:5*deg,j+1));
    end
    w23(:,j+1)=w23(:,j)-rho*(x(deg+1:2*deg,j+1)-x(2*deg+1:3*deg,j+1));
    w34(:,j+1)=w34(:,j)-rho*(x(2*deg+1:3*deg,j+1)-x(3*deg+1:4*deg,j+1));
    if(mal_node~=5)
        w45(:,j+1)=w45(:,j)-rho*(x(3*deg+1:4*deg,j+1)-x(4*deg+1:5*deg,j+1));
    end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %         Calculate E
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (FOUND==0) %if we have not found the malacious node
           E12=0.5*norm(x(deg+1:2*deg,j+1)-  x(1:deg,j+1));
           E15=0.5*norm(x(4*deg+1:5*deg,j+1) - x(1:deg,j+1)); % 5 is a successor of 1
           E23=0.5*norm(x(2*deg+1:3*deg,j+1) - x(deg+1:2*deg,j+1)); % 3 is a successor of 2
           E21=0.5*norm(x(1:deg,j+1) - x(deg+1:2*deg,j+1)); % 1 is a predecessor of 2
           E34=0.5*norm(x(3*deg+1:4*deg,j+1) -  x(2*deg+1:3*deg,j+1)); % 4 is a successor of 3
           E32=0.5*norm(x(deg+1:2*deg,j+1) -  x(2*deg+1:3*deg,j+1)); % 2 is a predecessor of 3
           E45=0.5*norm(x(4*deg+1:5*deg,j+1) -  x(3*deg+1:4*deg,j+1)); % 5 is a successor of 4
           E43=0.5*norm(x(2*deg+1:3*deg,j+1) -  x(3*deg+1:4*deg,j+1)); % 3 is a predecessor of 4
           E51=0.5*norm(x(1:deg,j+1) -  x(4*deg+1:5*deg,j+1)); % 1 is a predecessor of 5
           E54=0.5*norm(x(3*deg+1:4*deg,j+1) -  x(4*deg+1:5*deg,j+1)); % 4is a predecessor of 5

          %text=sprintf('E12=%d E23=%d E34=%d E51=%d E54=%d',E12,E23,E34,E51,E54);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %        Markov  Model Transition Function calculation
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          alpha=1/j;

           %Node 1, successors are 2,5
           %Here 2=area 2 , 3=area 5
           text=sprintf('E12=%d E15=%d',E12,E15);
           disp(text);
           E1(1,2)=DEFAULT+(alpha*E12)+(1-alpha)*E1(1,2);
           E1(1,3)=DEFAULT+(alpha*E15)+(1-alpha)*E1(1,3);
           E1(2,1)=E1(1,2);
           E1(2,3)=DEFAULT;
           E1(3,1)=E1(1,3);
           E1(3,2)=E1(2,3);
          
            display(E1);
            tot_sum=[0,0,0];
            tot_sum(1)=E1(1,1)+E1(1,2)+E1(1,3);
            tot_sum(2)=E1(2,1)+E1(2,2)+E1(2,3);
            tot_sum(3)=E1(3,1)+E1(3,2)+E1(3,3);
            display(tot_sum);
            for row=1:3
                for col=1:3
                  E1(row,col)=E1(row,col)/tot_sum(row);
                end
            end
          
          
          [ st,f,max_diff,error] = Stationary(E1,j,st_prev);
           st_prev=st;
           if(max_diff>global_diff)
               global_diff=max_diff;
           end
           if(error==1) 
                convergence_flag=TRUE;
                disp('1. convergence_flag=TRUE');
                break;
       	   end
%            if(max_diff>0.2)
%                disp('MALACIOUS NODE FOUND');
%                display(st);%global_diff(1)=max_diff;
%                  convergence_flag=TRUE;
%            end
           
       end%end of found
end  %  end %end of for j 
display(j);
% x_fault_free=PRECISION*x(:,ending_interation_num); 
% save('abc_fault_free.mat');
% x_b4_gaussian(:,j+1)=dummy_1\dummy_5;
 %            x((row-1)*deg+1:row*deg,j+1)

ending_interation_num_no_fault=j;
text=sprintf('Gasussain scale=0.01 \n %d %d \n', global_diff,j);
disp(text);
fprintf(fid, '%d %d \n', global_diff,j);
fclose(fid);

x_gaussian=PRECISION*x(:,ending_interation_num_no_fault); 

%display(x(:,ending_interation_num));
save('gaussian_0_001.mat');

