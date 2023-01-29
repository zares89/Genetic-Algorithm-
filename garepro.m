function [epot_tot_sel, struc_sel] = garepro(epot_totall,ptour,struc,M)

% Randomly select 4M structures

%rand_2M = randperm(2*M);
randnum = rand(2*M,1);


% Choosing parents
parents=zeros(2*M,1);
parfor i=1:2*M
    rand_2M = randperm(2*M,2);
    [esorted,isorted]= sort([epot_totall(rand_2M(1)) epot_totall(rand_2M(2))]);
    
    if randnum(i) >=ptour
        parents(i)= rand_2M(isorted(2));
    else
        parents(i)= rand_2M(isorted(1));
    end
end

% Condensing the structure/epot_tot matrices to the selected parents
struc_sel = struc(parents,:);
epot_tot_sel  = epot_totall(parents);
end