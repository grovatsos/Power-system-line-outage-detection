%Compute the bus incidence matrix of the network graph


function A = adjacency_matrix(line, nbus)
A=eye(nbus);
for i=1:size((line),1)
    A(line(i,1),line(i,2))=1;
    A(line(i,2),line(i,1))=1;
end

