 function  struc= finestruc( struc,n,n_new,natom)
 
 

n_struc = size(struc,1);
struc_new = [];
%epot_tot = zeros(n_struc,1);
% f_cost = 4*eps*((sig/dist)^12-(sig/dist)^6);
parfor i=1:n_struc
    tmp_struc = [];
    for j=1:natom
     %   [i j]
    %    for k=j+1:natom
            % Binary to decimal
            coordfx = (1/2^n*bi2de(struc(i,(j-1)*3*n+1:(j-1)*3*n+n)));
            
            coordfy = (1/2^n*bi2de(struc(i,(j-1)*3*n+n+1:(j-1)*3*n+2*n)));
            coordfz = (1/2^n*bi2de(struc(i,(j-1)*3*n+2*n+1:(j-1)*3*n+3*n)));
            coordbx = de2bi(mod(round(coordfx*2^(n_new)),2^(n_new)),n_new);
            coordby = de2bi(mod(round(coordfy*2^(n_new)),2^(n_new)),n_new);
            coordbz = de2bi(mod(round(coordfz*2^(n_new)),2^(n_new)),n_new);
            tmp_struc = [tmp_struc coordbx coordby coordbz];
          %  coordkx = bi2de(struc(i,(k-1)*3*n+1:(k-1)*3*n+n));
          %  coordky = bi2de(struc(i,(k-1)*3*n+n+1:(k-1)*3*n+2*n));
          %  coordkz = bi2de(struc(i,(k-1)*3*n+2*n+1:(k-1)*3*n+3*n));
         %   % Return the distance components to fractional coordinate
         %   fracdx   = 1/(2^n)* (coordjx-coordkx);
        %    fracdy   = 1/(2^n)* (coordjy-coordky);
        %    fracdz   = 1/(2^n)* (coordjz-coordkz);
            % Lowest periodic distance
        %    fracdx   = fracdx - round(fracdx);
       %     fracdy   = fracdy - round(fracdy);
        %    fracdz   = fracdz - round(fracdz);
            
        %    distfrac = sqrt(fracdx^2+fracdy^2+fracdz^2);
            
            % Poitential energy calculation
        %    epot = 4*eps*((sigma/(cella * distfrac))^12-(sigma/(cella * distfrac))^6);
        %    epot_tot(i)= epot_tot(i) + epot;
   %     end
    end
    struc_new = [struc_new;tmp_struc];
end
struc = struc_new;
end