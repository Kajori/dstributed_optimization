% 3 phase fault  -no
% Gaussian Noise - YES
% malacious node - if present node 5
% every node will receive values from evry other node
% there will be a central guy.. trim the top f and bottom f values 
% The value added by the adversaty is epsilon less then the values proposed
% by the other guys.. hence his value never gets thrown away
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
N=100;
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
x_5(:,:,1)=ones(NO_AREA-1,deg,1); %it is the value sent by node 5 to different nodes
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
    
    outliers=zeros(NO_AREA,2);
    for row=1:NO_AREA
       
        [left,right]=cal_left_right(row);
        left_val= cal_outlier_individual(x,x_5,j,row,NO_AREA,deg,left);
        right_val= cal_outlier_individual(x,x_5,j,row,NO_AREA,deg,right);
        Pre_succ=2;
        
        dummy=reshape(H_pdc(row,1:TT,:),TT,deg);
        dummy_0=dummy.';
        dummy_01=dummy_0*dummy;
        dummy_1=dummy_01+rho*Pre_succ*eye(size(dummy,2));
        dummy_3=reshape(C_pdc(row,1:TT,:),TT,1);
        dummy_4=dummy'*dummy_3;
       
        if(row==1)
            dummy_5=dummy_4+(w12_1(:,j) + w15_1(:,j))+rho*(left_val+right_val);
            x((row-1)*deg+1:row*deg,j+1) =dummy_1\dummy_5;
         elseif(row==2)
              dummy_5=dummy_4+(w23_2(:,j) - w12_2(:,j))+rho*(left_val+right_val);
              x((row-1)*deg+1:row*deg,j+1) =dummy_1\dummy_5;
              w12_2(:,j+1)=w12_2(:,j)-rho*(left_val-x(deg+1:2*deg,j+1));   % w12 update as 2 is the successor of 1
        elseif(row==3)
             dummy_5=dummy_4+(w34_3(:,j) - w23_3(:,j))+rho*(left_val+right_val);
             x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
             w23_3(:,j+1)=w23_3(:,j)-rho*(left_val-x(2*deg+1:3*deg,j+1));   % w23 update as 2 is the predecessor of 3
       elseif(row==4)
            dummy_5=dummy_4+(w45_4(:,j) - w34_4(:,j))+rho*(left_val+right_val);
            x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
            w34_4(:,j+1)=w34_4(:,j)-rho*(left_val-x(3*deg+1:4*deg,j+1));   
       elseif ( row==5)
            dummy_5=dummy_4+((-1)*w15_5(:,j) - w45_5(:,j))+rho*(left_val+right_val);
            x((row-1)*deg+1:row*deg,j+1)=dummy_1\dummy_5;
            x_5=clever_adversary(row,x_5,NO_AREA,x,j+1,deg);
            %x_5=good_adversary(row,x_5,NO_AREA,x,j+1,deg);
            w15_5(:,j+1)=w15_5(:,j)-rho*(left_val-x(4*deg+1:5*deg,j+1));        % w14 update as 1 is the predecessor of 4
            w45_5(:,j+1)=w45_5(:,j)-rho*(right_val-x(4*deg+1:5*deg,j+1));        % w14 update as 1 is the predecessor of 4
       end %end of if
   end %end of for row
   
    convergence_flag=TRUE;
    for row=1:NO_AREA*deg
        if ((abs(x(row,j+1)-x(row,j)))>0.000001)
          if(j>1000) disp(sprintf('x(%d,%d)-x(%d,%d)=' row,j+1,row,j)); convergence_flag=FALSE; end
        end
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %         after receiving from all successor
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     row=1;
     [left,right]=cal_left_right(row);
     left_val= cal_outlier_individual(x,x_5,j+1,row,NO_AREA,deg,left);
     right_val= cal_outlier_individual(x,x_5,j+1,row,NO_AREA,deg,right);
     w15_1(:,j+1)=w15_1(:,j)-rho*(x(1:deg,j+1)-left_val);
     w12_1(:,j+1)=w12_1(:,j)-rho*(x(1:deg,j+1)-right_val);
     
     row=2;
     [~,right]=cal_left_right(row);
     right_val= cal_outlier_individual(x,x_5,j+1,row,NO_AREA,deg,right);
     w23_2(:,j+1)=w23_2(:,j)-rho*(x(deg+1:2*deg,j+1)-right_val);
     
     row=3;
     [~,right]=cal_left_right(row);
     right_val= cal_outlier_individual(x,x_5,j+1,row,NO_AREA,deg,right);
     w34_3(:,j+1)=w34_3(:,j)-rho*(x(2*deg+1:3*deg,j+1)-right_val);
     
     row=4;
     [~,right]=cal_left_right(row);
     right_val= cal_outlier_individual(x,x_5,j+1,row,NO_AREA,deg,right);
     w45_4(:,j+1)=w45_4(:,j)-rho*(x(3*deg+1:4*deg,j+1)-right_val);
end  %  end %end of for j 
display(j);
disp('row == 1');
disp(x(1:deg,j));
disp('row == 2');
disp(x(deg+1:2*deg,j));
%display(x(:,ending_interation_num));
save('fault_no_detection_3per.mat');

