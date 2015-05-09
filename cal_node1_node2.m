% there if the neighbours are malacious its ets the node 1 and node 2
% respectively fot the left and right neighbours


function [n1,n2] = cal_node1_node2(row,mal_node_1,mal_node_2)

    n1=0;
    n2=0;
    if(row==1 && mal_node_1==5)       n1=1;
    elseif(row==2 && mal_node_1==1)   n1=1;
    elseif(row==3 && mal_node_1==2)   n1=1;
    elseif(row==4 && mal_node_1==3)   n1=1;
    elseif(row==5 && mal_node_1==4)   n1=1;
    end
    
    if(row==1 && mal_node_2==2)       n2=1;
    elseif(row==2 && mal_node_2==3)   n2=1;
    elseif(row==3 && mal_node_2==4)   n2=1;
    elseif(row==4 && mal_node_2==5)   n2=1;
    elseif(row==5 && mal_node_2==1)   n2=1;
    end
    
    
