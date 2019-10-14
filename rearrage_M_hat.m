function [M_hat_rearranged] = rearrage_M_hat(M_hat,gen_indices,load_indices,slack_index)
%rearrage the rows of M0_hat so that the buses for theta corresponds to original bus indices
temp=M_hat*0;
for i=1:length(M_hat)
        temp_row=M_hat(i,:)   
        i_new=orig_bus(i+1, gen_indices, load_indices, slack_index);
        temp(i_new-1,:)=temp_row;    
end
M_hat_rearranged=temp;
end

