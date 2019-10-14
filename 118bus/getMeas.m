function A = getMeas(mpc,Psig,line)

    define_constants;
    baseMVA = mpc.baseMVA;
    bus = mpc.bus;
    nbus = size(bus,1);
    type = round(bus(:,2));
    
    P = bus(:,3) + Psig*baseMVA*randn(nbus,1);
    mpc.branch = line;
    mpc.bus(:,3) = P;
    
    mpopt = mpoption('PF_ALG', 1, 'VERBOSE', 0, 'OUT_ALL', 0);
    results = runpf(mpc, mpopt);
    A = (results.bus(type~=3,VA) - results.bus(type==3,VA))*pi/180;

end