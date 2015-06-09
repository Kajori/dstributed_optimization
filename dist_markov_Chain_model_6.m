%3 phase fault  -no
%Gaussian Noise - YES
% malacious node - if present node 5
% there will be a central guy.. trim the top f and bottom f values 
%here f =1
MAX_DEV = 1e-6;

clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Defining the constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NO_AREA=5;
PMU_PER_AREA=3;
SIZE_OF_voltage_data=1000;
NO_OF_EIGEN_VALUES=20;

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
%display('Gaussian Noise no 3 phase fault');

rho=10^-9;
disp('-----------------------------------------------');
tic;

TT=PMU_PER_AREA*m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Loading data from the input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
for row=1:NO_AREA
     for col=1:PMU_PER_AREA
           file_name=sprintf('no_3_phase_fault_area_%d_line_%d.txt',row,area_line_no(row,col));
           Y_org(row,col,:)=importdata(file_name,'\n');
      end %end of for col
end %end of for row
%disp(sprintf( ' size(Y_org) =  %d ',size(Y_org)));
disp(' Loaded data from the input files ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Write statistics to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('stat.txt', 'a+');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %         adding Gausian Noise  to the data
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCALE = 0.001; %0.01 % noise
for row=1:NO_AREA 
    for col=1:PMU_PER_AREA
             noise = randn(1); % noise with mean=0 and std=1;
             Y(row,col,:)= ( 1+  noise* SCALE )* Y_org(row,col,:);
    end %end of for col
end %end of for row        
%disp(sprintf( ' size(Y) = %d',size(Y)));

     
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
     
     %disp(text);
     %disp(size(H_pdc));
     %disp(size(C_pdc));
 end %end of for row



for rn=0:NO_AREA %rn means the removed node
w=zeros(NO_AREA*deg,1);
x(:,1)=ones(NO_AREA*deg,1); %Everyone has equal initial x
z(:,1)=zeros(deg,1);
for j=1:10000
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         if we have not reached convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text=sprintf(' -------------- %d --------------',j);
    disp( text);
    for row=1:NO_AREA
        if(row~=rn)
            H_r=reshape(H_pdc(row,1:TT,:),TT,deg);
            C_r=reshape(C_pdc(row,1:TT,:),TT,1);
            temp=(H_r'*C_r)-w((row-1)*deg+1:row*deg,j)+rho*z(:,j);
            x((row-1)*deg+1:row*deg,j+1)=(H_r'*H_r+rho*eye(size(H_r,2)))\temp;
            if(row==5) x((row-1)*deg+1:row*deg,j+1)=1.05*x((row-1)*deg+1:row*deg,j+1);
            end
        end
   end %end row
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Remove one node from the calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   temp=zeros(deg,1);
   count=0;
   for row=1:NO_AREA
       if(row~=rn)
            temp=temp+x((row-1)*deg+1:row*deg,j+1);
            count=count+1;
       end
   end
   temp=temp/count;
   z(:,j+1)=transpose(temp');
   
   for row=1:NO_AREA
       if(row~=rn)
           w((row-1)*deg+1:row*deg,j+1)= w((row-1)*deg+1:row*deg,j)+rho*( x((row-1)*deg+1:row*deg,j+1)-z(:,j+1)); 
       end
   end %end row
    
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %         Check if we have reached convergence
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   convergence_flag=TRUE;
   for row=1:NO_AREA*deg
        if(row~=rn)
            if ((abs(x(row,j+1)-x(row,j)))>0.0000001)
                convergence_flag=FALSE;
            end
        end
   end
   if( convergence_flag==TRUE)
         display(sprintf('Convergence Reached at %d ',j));
         break;
    end 
end %end of j

if (rn==0) x_0=x(:,j);
elseif (rn==1) x_1=x(:,j);
elseif (rn==2) x_2=x(:,j);
elseif (rn==3) x_3=x(:,j);
elseif (rn==4) x_4=x(:,j);
elseif (rn==5) x_5=x(:,j);
end    

end %for rn=1:NO_AREA

disp(' \n ');
disp('%%%%%%%%%%%%%% for area 1 %%%%%%%%%%%%%%%%');
for row=1:deg
    disp(sprintf(' %s %s  %s  %s  %s',x_0(row),x_2(row),x_3(row),x_4(row),x_5(row)));
end

disp(' \n ');
disp(' %%%%%%%%%%%%%% for area 2 %%%%%%%%%%%%%%%%');
for row=deg+1:2*deg
    disp(sprintf(' %s %s  %s  %s  %s',x_0(row),x_1(row),x_3(row),x_4(row),x_5(row)));
end

disp(' \n ');
disp('%%%%%%%%%%%%%% for area 3 %%%%%%%%%%%%%%%%');
for row=2*deg+1:3*deg
    disp(sprintf(' %s  %s  %s  %s',x_0(row),x_1(row),x_2(row),x_4(row),x_5(row)));
end

disp(' \n ');
disp(' %%%%%%%%%%%%%% for area 4 %%%%%%%%%%%%%%%%');
for row=3*deg+1:4*deg
    disp(sprintf(' %s %s  %s  %s  %s',x_0(row),x_1(row),x_2(row),x_3(row),x_5(row)));
end

disp(' \n ');
disp(' %%%%%%%%%%%%%% for area 5 %%%%%%%%%%%%%%%%');
for row=4*deg+1:5*deg
    disp(sprintf(' %s %s  %s  %s  %s',x_0(row),x_1(row),x_2(row),x_3(row),x_4(row)));
end

