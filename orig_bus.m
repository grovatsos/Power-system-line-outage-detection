function [bus_orig] = orig_bus(bus_dsb, gen_indices, load_indices, slack_index)
% Input: Distributed slack bus values, Output: Original buses

    if (bus_dsb==1)
         bus_orig=slack_index;
    elseif (bus_dsb<=numel(gen_indices)+1)
        bus_orig=gen_indices(bus_dsb-1);
    else
         bus_orig=load_indices(bus_dsb-numel(gen_indices)-1);
    end
end

