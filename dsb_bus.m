function [bus_dsb] = dsb_bus(bus_orig, gen_indices, load_indices, slack_index)
% Input: Old bus values, Output: New bus values for distributed slack bus power flows
    s=find(slack_index==bus_orig);
    g=find(gen_indices==bus_orig);
    l=find(load_indices==bus_orig);
        if s==1
            bus_dsb=1;
        elseif (g>0)
            bus_dsb=g+1;%bus for new rearranged system where generator buses are stacked first then load buses +1 for slack
        else
            bus_dsb=numel(gen_indices)+l+1;
        end
end

