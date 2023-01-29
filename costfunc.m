 function [epot_tot,bestn,bestepot] = costfunc( natom, struc, sigma, eps, cella, n,N)

%%Evaluate cost function
n_struc = size(struc,1);
epot_tot = zeros(n_struc,1);
% f_cost = 4*eps*((sig/dist)^12-(sig/dist)^6);
parfor i=1:n_struc
    for j=1:natom
        for k=j+1:natom
            % Binary to decimal
            coordjx = bi2de(struc(i,(j-1)*3*n+1:(j-1)*3*n+n));
            
            coordjy = bi2de(struc(i,(j-1)*3*n+n+1:(j-1)*3*n+2*n));
            coordjz = bi2de(struc(i,(j-1)*3*n+2*n+1:(j-1)*3*n+3*n));
            coordkx = bi2de(struc(i,(k-1)*3*n+1:(k-1)*3*n+n));
            coordky = bi2de(struc(i,(k-1)*3*n+n+1:(k-1)*3*n+2*n));
            coordkz = bi2de(struc(i,(k-1)*3*n+2*n+1:(k-1)*3*n+3*n));
            % Return the distance components to fractional coordinate
            fracdx   = 1/(2^n)* (coordjx-coordkx);
            fracdy   = 1/(2^n)* (coordjy-coordky);
            fracdz   = 1/(2^n)* (coordjz-coordkz);
            % Lowest periodic distance- Uncomment for PBC consideration
            
          %  fracdx   = fracdx - round(fracdx);
          %  fracdy   = fracdy - round(fracdy);
          %  fracdz   = fracdz - round(fracdz);
            
            distfrac = sqrt(fracdx^2+fracdy^2+fracdz^2);
            
            % Poitential energy calculation
            epot = 4*eps*((sigma/(cella * distfrac))^12-(sigma/(cella * distfrac))^6);
            epot_tot(i)= epot_tot(i) + epot;
        end
    end
end
[pot_sorted, id_sorted] = sort(epot_tot);
bestn = struc(id_sorted(1:N),:);
bestepot = pot_sorted(1:N);
 end