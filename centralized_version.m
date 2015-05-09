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

text = sprintf(' size(H_pmu) size(C_pmu) ');
%disp(text);
%disp(size(H_pmu));
%disp(size(C_pmu));

H_pdc=zeros(NO_AREA,PMU_PER_AREA*m,deg);
C_pdc=zeros(NO_AREA,PMU_PER_AREA*m);
%display('Gaussian Noise no 3 phase fault');

rho=10^-9;
N=10000;
disp('-----------------------------------------------');
tic;

TT=PMU_PER_AREA*m;
SCALE = 0.001; %0.01 % noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Loading data from the input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
for row=1:NO_AREA
     for col=1:PMU_PER_AREA
           file_name=sprintf('no_3_phase_fault_area_%d_line_%d.txt',row,area_line_no(row,col));
           Y_org(row,col,:)=importdata(file_name,'\n');
      end %end of for col
end %end of for row
text = sprintf( ' size(Y_org) =  ');
%disp(text);
%disp(size(Y_org));
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
text = sprintf( ' size(Y) =');
%disp(text);
%disp(size(Y));
     
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
     
     text = sprintf( '1. row = %d size(H_pdc) size(C_pdc)=%d ',row);
     %disp(text);
     %disp(size(H_pdc));
     %disp(size(C_pdc));
 end %end of for row

w=zeros(NO_AREA*deg,1);
x(:,1)=ones(NO_AREA*deg,1); %Everyone has equal initial x
z(:,1)=zeros(deg,1);
sum_violation=zeros(1,NO_AREA);
convergence_flag=FALSE;
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
    for row=1:NO_AREA 
        H_r=reshape(H_pdc(row,1:TT,:),TT,deg);
        C_r=reshape(C_pdc(row,1:TT,:),TT,1);
        temp=(H_r'*C_r)-w((row-1)*deg+1:row*deg,j)+rho*z(:,j);
        x((row-1)*deg+1:row*deg,j+1)=(H_r'*H_r+rho*eye(size(H_r,2)))\temp;
        %disp ( ' ******************* ');
        %disp(row);
%         disp(size(temp));
%         disp(size(temp_2));
%         disp(size(temp_3));
%         disp(size(temp_4));
   end %end row
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         find the top f and bottom f outliers
    %         max_sum and min_sum are the top and bottom outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for row=1:NO_AREA
       sum_temp=0;
       for col=1:NO_AREA
            sum_temp=sum_temp + 0.5*norm(x((row-1)*deg+1:row*deg,j+1)-  x((col-1)*deg+1:col*deg,j+1));
       end
       sum_violation(1,row)=sum_temp;
   end
   
   [max_sum, max_index ]=max(sum_violation);
   [min_sum, min_index ]=min(sum_violation);
   
   temp=zeros(deg,1);
   for row=1:NO_AREA
       if(row~=max_index && row~=min_index)
            temp=temp+x((row-1)*deg+1:row*deg,j+1);
       end
    end
    z(:,j+1)=transpose(temp');
    %disp(size(z));
    disp('sum_violation =');
    disp(sum_violation);
    text = sprintf (' max_index = %d min_index = %d ',max_index,min_index);
    disp(text);
   for row=1:NO_AREA
       if(row~=max_index && row~=min_index)
            w((row-1)*deg+1:row*deg,j+1)= w((row-1)*deg+1:row*deg,j)+rho*( x((row-1)*deg+1:row*deg,j+1)-z(:,j+1)); 
       else 
            w((row-1)*deg+1:row*deg,j+1)= w((row-1)*deg+1:row*deg,j);
       end
   end %end row
    
   
%    disp('max_sum');
%    disp(max_sum);
%    disp('index');
%    disp(index);
%    disp('sum_violation');
%    disp(sum_violation);
    
       
            
           
   
   
   convergence_flag=TRUE;
   for row=1:NO_AREA*deg
        if ((abs(x(row,j+1)-x(row,j)))>0.000001)
          convergence_flag=FALSE;
        end
   end
    
end
%print x(:,j)
