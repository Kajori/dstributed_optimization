% 3 phase fault  -no
% Gaussian Noise - YES
% malacious node - if present node 5
% every node will receive values from evry other node
% there will be a central guy.. trim the top f and bottom f values 
% here f =1
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
SCALE = 0.001; %0.1 % noise

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

w=zeros(deg,1); 
w12_1(:,1)=w; w15_1(:,1)=w; 
w23_2(:,1)=w; w34_3(:,1)=w; w45_4(:,1)=w; %Everyone has equal initial w_ij
w12_2(:,1)=w; w15_5(:,1)=w; 
w23_3(:,1)=w; w34_4(:,1)=w; w45_5(:,1)=w; %Everyone has equal initial w_ij
x(:,1)=ones(NO_AREA*deg,1);

convergence_flag=FALSE;


for j=1:N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Check if we have reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( convergence_flag==TRUE)
        display(sprintf('Convergence Reached at %d ',j));
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         if we have not reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(sprintf(' -------------- %d --------------',j));
    
    for row=1:NO_AREA
       
        
        [ mal_node_1,mal_node_2,Pre_succ]= cal_outlier(x,j,row,NO_AREA,deg);
        mal_node_1=-1;
        mal_node_2=-1;
        Pre_succ=2;
        
        dummy=reshape(H_pdc(row,1:TT,:),TT,deg);
        dummy_0=dummy.';
        dummy_01=dummy_0*dummy;
        dummy_1=dummy_01+rho*Pre_succ*eye(size(dummy,2));
        dummy_3=reshape(C_pdc(row,1:TT,:),TT,1);
        dummy_4=dummy'*dummy_3;
       
        
       
    
        if(row==1)
            if  (mal_node_1==5 && mal_node_2==2) 
                x((row-1)*deg+1:row*deg,j+1)= x((row-1)*deg+1:row*deg,j);
            else
                 if  (mal_node_1~=5 && mal_node_2~=2)
                    dummy_5=dummy_4+(w12_1(:,j) + w15_1(:,j))+rho*(x(deg+1:2*deg,j)+x(4*deg+1:5*deg,j));
                elseif  (mal_node_1~=5 && mal_node_2==2)
                    dummy_5=dummy_4+(w15_1(:,j))+rho*(x(4*deg+1:5*deg,j));
                elseif  (mal_node_1==5 && mal_node_2~=2)
                    dummy_5=dummy_4+(w12_1(:,j))+rho*(x(deg+1:2*deg,j));
                 end
                x((row-1)*deg+1:row*deg,j+1) =dummy_1\dummy_5;
            end
            
        elseif(row==2)
            if  (mal_node_1==1 && mal_node_2==3) x((row-1)*deg+1:row*deg,j+1)=x((row-1)*deg+1:row*deg,j);
            else    
                if  (mal_node_1~=1 && mal_node_2~=3)
                    dummy_5=dummy_4+(w23_2(:,j) - w12_2(:,j))+rho*(x(1:deg,j+1)+x(2*deg+1:3*deg,j));
                elseif  (mal_node_1~=1 && mal_node_2==3)
                    dummy_5=dummy_4 - w12_2(:,j)+rho*(x(1:deg,j+1));
                elseif  (mal_node_1==1 && mal_node_2~=3)
                    dummy_5=dummy_4+(w23_2(:,j))+rho*(x(2*deg+1:3*deg,j));
                end 
                x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
            end  
            
            if  (mal_node_1~=1)
                w12_2(:,j+1)=w12_2(:,j)-rho*(x(1:deg,j+1)-x(deg+1:2*deg,j+1));   % w12 update as 2 is the successor of 1
            else
                w12_2(:,j+1)=w12_2(:,j);
            end
            
            
        elseif(row==3)
             if  (mal_node_1==2 && mal_node_2==4) 
                x((row-1)*deg+1:row*deg,j+1)= x((row-1)*deg+1:row*deg,j)
             else
                if  (mal_node_1~=2 && mal_node_2~=4)
                    dummy_5=dummy_4+(w34_3(:,j) - w23_3(:,j))+rho*(x(deg+1:2*deg,j+1)+x(3*deg+1:4*deg,j));
                elseif  (mal_node_1~=2 && mal_node_2==4)
                    dummy_5=dummy_4+( - w23_3(:,j))+rho*(x(deg+1:2*deg,j+1));
                elseif  (mal_node_1==2 && mal_node_2~=4)
                    dummy_5=dummy_4+(w34_3(:,j) )+rho*(x(3*deg+1:4*deg,j));
                end
                x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
             end
             
             if  (mal_node_1==2)
                 w23_3(:,j+1)=w23_3(:,j);
             else
                w23_3(:,j+1)=w23_3(:,j)-rho*(x(deg+1:2*deg,j+1)-x(2*deg+1:3*deg,j+1));   % w23 update as 2 is the predecessor of 3
             end
             
        elseif(row==4)
            if  (mal_node_1==3 && mal_node_2==5) x((row-1)*deg+1:row*deg,j+1)=x((row-1)*deg+1:row*deg,j);
            else
                if  (mal_node_1~=3 && mal_node_2~=5)
                    dummy_5=dummy_4+(w45_4(:,j) - w34_4(:,j))+rho*(x(2*deg+1:3*deg,j+1)+x(4*deg+1:5*deg,j));
                elseif  (mal_node_1~=3 && mal_node_2==5)
                    dummy_5=dummy_4 - w34_4(:,j)+rho*(x(2*deg+1:3*deg,j+1));
                elseif  (mal_node_1==3 && mal_node_2~=5)
                    dummy_5=dummy_4+(w45_4(:,j))+rho*(x(4*deg+1:5*deg,j));
                end
                x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
            end
            if  (mal_node_1~=3) % w34 update as 3 is the predecessor of 4
                w34_4(:,j+1)=w34_4(:,j)-rho*(x(2*deg+1:3*deg,j+1)-x(3*deg+1:4*deg,j+1));   
            else
                w34_4(:,j+1)=w34_4(:,j);
            end
                
            
      elseif ( row==5)
           if  (mal_node_1==4 && mal_node_2==1) x((row-1)*deg+1:row*deg,j+1)=x((row-1)*deg+1:row*deg,j);
           else
               if (mal_node_1~=4 && mal_node_2~=1)
                     dummy_5=dummy_4+(-w15_5(:,j) - w45_5(:,j))+rho*(x(1:deg,j+1)+x(3*deg+1:4*deg,j+1));
               elseif (mal_node_1~=4 && mal_node_2==1)
                     dummy_5=dummy_4+(- w45_5(:,j))+rho*(x(3*deg+1:4*deg,j+1));
               elseif (mal_node_1==4 && mal_node_2~=1)
                     dummy_5=dummy_4+(-w15_5(:,j))+rho*(x(1:deg,j+1));
               end   
              x((row-1)*deg+1:row*deg,j+1)=1.03*dummy_1\dummy_5;
           end
           
           if ( mal_node_2~=1)
                w15_5(:,j+1)=w15_5(:,j)-rho*(x(1:deg,j+1)-x(4*deg+1:5*deg,j+1));        % w14 update as 1 is the predecessor of 4
           else
                 w15_5(:,j+1)=w15_5(:,j);
           end
           if (mal_node_1~=4)
                w45_5(:,j+1)=w45_5(:,j)-rho*(x(3*deg+1:4*deg,j+1)-x(4*deg+1:5*deg,j+1));        % w14 update as 1 is the predecessor of 4
           else
                w45_5(:,j+1)=w45_5(:,j);
           end
        end %end of if
   end %end of for row
   
   convergence_flag=TRUE;
    for row=1:NO_AREA*deg
        if ((abs(x(row,j+1)-x(row,j)))>0.000001)
          convergence_flag=FALSE;
        end
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %         after receiving from all successor
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w12_1(:,j+1)=w12_1(:,j)-rho*(x(1:deg,j+1)-x(deg+1:2*deg,j+1));
    if(mal_node~=5)
          w15_1(:,j+1)=w15_1(:,j)-rho*(x(1:deg,j+1)-x(4*deg+1:5*deg,j+1));
    end
    w23_2(:,j+1)=w23_2(:,j)-rho*(x(deg+1:2*deg,j+1)-x(2*deg+1:3*deg,j+1));
    w34_3(:,j+1)=w34_3(:,j)-rho*(x(2*deg+1:3*deg,j+1)-x(3*deg+1:4*deg,j+1));
    if(mal_node~=5)
        w45_4(:,j+1)=w45_4(:,j)-rho*(x(3*deg+1:4*deg,j+1)-x(4*deg+1:5*deg,j+1));
    end
end  %  end %end of for j 
display(j);

x_fault_no_detection=x;
%display(x(:,ending_interation_num));
save('fault_no_detection_3per.mat');

