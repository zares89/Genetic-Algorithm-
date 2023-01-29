function genetic_material

fileenergy = 'energy-9-cell-9-pbc.gin';
natom = 9; 
n_init =1000; % no of initial structures
cella = 12; % cell dimension in Angstrom
nd_init = 4;   % Number of points on one dimension is 2^n+1
M = 50;
N = 2;L = 2; Nf = 5;%number of best candidates- number of foreign candidates- Number of final structures for optimization
max_iter = 500; %(number of iteration)
%max_ngrid =12; 
% MAterial properties (LJ12-6 parameters, sigma should be dimensionless)
eps = 300; sigma =3.5;
% Genetic parameters
ptour = 0.9;pcross = 0.5;% pmute = 0.0001;
%% Define the grid
%gridx = linspace(0,1,2^n+1)*cella;
%gridy = gridx;gridz = gridx;
% fractional grid
% gridfx = linspace(0,1,2^n+1);
% gridfy = gridfx;gridfz = gridfx;
% % Convert grid lines to binary
% gridfx(end)=[];gridfy(end)=[];gridfz(end)=[];
% gridbx = de2bi( gridfx*2^n );gridby = gridbx; gridbz = gridbx;
% % assign binary id to each point on the grid
% lgx = length(gridfx); lgy = length(gridfy) ; lgz = length(gridfz);
% pointid = zeros(lgx*lgy*lgz,n*3);

% COncentrating all to make the DNA
% cnt=0;
% for i=1:lgx
%     for j=1:lgy
%         for k=1:lgz
%             cnt= cnt+1;
%             pointid(cnt,:) = [gridbx(i,:) gridby(j,:) gridbz(k,:)];
%            
%         end
%     end
% end


%% Construction of initial random structures 
n =nd_init;
[struc_iter,struc_init] = initstruc(natom, n, n_init, M);
%% Iterating 
Mnew = M;

%struc_iter = struc;
pot_array=[];
tol_grid= 1e-4;
pmute = 1/ (3*natom*n);
for n_iter=1:max_iter
    % Call costfunc to evaluate the cost function
    n_iter
     if n_iter > 2 %%&& n <  max_ngrid
        if  mod(n_iter,20)==0 %%&& abs(pot_array(end,1)-pot_array(end-1,1)) < tol_grid
          %  epot_tot_test = costfunc( natom, struc_iter, sigma, eps, cella, n,N);
         n_new = n+1;
          struc_iter= finestruc( struc_iter,n,n_new,natom);
            struc_init= finestruc( struc_init,n,n_new,natom);
            n=n+1;
         %   epot_tot_test2 = costfunc( natom, struc_iter, sigma, eps, cella, n,N);
            pmute = 1/ (3*natom*n);
        end
    % elseif n== max_ngrid
    %     n_new = nd_init;
     %       struc_iter= finestruc( struc_iter,n,n_new,natom);
     %       struc_init= finestruc( struc_init,n,n_new,natom);
     %    n = n_new;
    end
[epot_totall,bestn,bestepot] = costfunc( natom, struc_iter, sigma, eps, cella, n,N);
    % Call reproduction function to pair parents
    [epot_tot_sel, struc_sel] = garepro(epot_totall,ptour,struc_iter,Mnew);
    % Crossing
    for i=1:Mnew
        randpar = rand(1);
        if randpar < pcross
            tmppar= struc_sel((i-1)*2+1,:);
            randlength = randi([1 natom*3*n],1);
            randpos1    = randi([1 natom*3*n-randlength+1],1);
            randpos2    = randi([1 natom*3*n-randlength+1],1);
            struc_sel((i-1)*2+1,randpos1:randpos1+randlength-1)=struc_sel((i-1)*2+2,randpos2:randpos2+randlength-1);
            struc_sel((i-1)*2+2,randpos2:randpos2+randlength-1)=tmppar(randpos1:randpos1+randlength-1);
            
        end
    end
    % Mutation
    for i=1:Mnew
        randpar = rand(1);
        if randpar < pmute
            % tmppar= struc(parents((i-1)*2+1),:);
            % Have no idea what is the best value for max number of swaps
            
            % Mutation in parent 1
            zerosz    = sum(struc_sel((i-1)*2+1,:)==0);
            onesz     = 3*n*natom- zerosz;
            randnswap0 = randi([1 zerosz],1);
            randnswap1 = randi([1  onesz],1);
            [zz,zerosid]    = find(struc_sel((i-1)*2+1,:)==0);
            [oo,onesid]    = find(struc_sel((i-1)*2+1,:)==1);
            randzeros = randi([1 zerosz],randnswap0,1);
            randones   = randi([1 onesz],randnswap1,1);
            struc_sel((i-1)*2+1,zerosid(randzeros)) =  1;
            struc_sel((i-1)*2+1,onesid(randones)) =  0;
            % Mutation in parent 2
           
            zerosz    = sum(struc_sel((i-1)*2+2,:)==0);
            onesz     = 3*n*natom- zerosz;
             randnswap0 = randi([1 zerosz],1);
            randnswap1 = randi([1 onesz],1);
            [zz,zerosid]    = find(struc_sel((i-1)*2+2,:)==0);
            [oo,onesid]    = find(struc_sel((i-1)*2+2,:)==1);
            randzeros = randi([1 zerosz],randnswap0,1);
            randones   = randi([1 onesz],randnswap1,1);
            struc_sel((i-1)*2+2,zerosid(randzeros)) =  1;
            struc_sel((i-1)*2+2,onesid(randones)) =  0;
            
        end
    end
    
    %% Block for calculating energy of children
    
    [epot_children,bestdummy,bestedummy] = costfunc( natom, struc_sel, sigma, eps, cella, n,N);

            %%
    bestedummy'
    pot_array=[pot_array;bestedummy'];
  %  pot_array(end,:)
    Lnew   = randperm(size(struc_init,1),L);
    struc_L = struc_init(Lnew,:);
    struc_iter = [struc_sel;bestn;struc_L];
    
    % Shuffle all the new parents and calculate new cost functions
    struc_iter = struc_iter(randperm(2*Mnew+N+L),:); % Check the size here
    Mnew = size(struc_iter,1)*0.5;
    
   
end
struc_final = struc_sel;
[epot_final,bestnfinal,bestepotfinal] = costfunc( natom, struc_final, sigma, eps, cella, n,Nf);


% FOrcefield transformation
A = 4*eps*sigma^12;B = 4* eps * sigma^6;
for i=1:Nf
    coordx = [];
    coordy = [];
    coordz = [];
    for j=1:natom
        coordx = [coordx;1/(2^n)*bi2de(bestnfinal(i,(j-1)*3*n+1:(j-1)*3*n+n))];
        coordy = [coordy;1/(2^n)*bi2de(bestnfinal(i,(j-1)*3*n+n+1:(j-1)*3*n+2*n))];
        coordz = [coordz;1/(2^n)*bi2de(bestnfinal(i,(j-1)*3*n+2*n+1:(j-1)*3*n+3*n))];
    end
    %% Write to Gulp runfile
         filename = ['struc-best-' num2str(i) '.gin'];
         gulpfile = fopen(filename,'w');
        fprintf(gulpfile,'opti conv noelectro\n');
   %     fprintf(gulpfile,'cell\r\n');
   %     fprintf(gulpfile,'%4.8f %4.8f %4.8f 90 90 90\r\n',cella,cella,cella);
        fprintf(gulpfile,'cart\r\n');
        for k=1:length(coordx)
            
        fprintf(gulpfile,'Al core %4.8f %4.8f %4.8f\r\n',cella*coordx(k),cella*coordy(k),cella*coordz(k));
        end
        fprintf(gulpfile,'lennard 12  6\n');
        fprintf(gulpfile,'Al    core Al    core %4.8f  %4.8f      0.000 10.0\n',A,B);
        dumptxt =['dump every 10000 best-' num2str(i) '.grs cart \r\n'];
        xyztxt =['output movie 1000 xyz bestga-' num2str(i) '.xyz \r\n'];
        fprintf(gulpfile,dumptxt);
        fprintf(gulpfile, xyztxt);
        
end
fclose(gulpfile);
%% Write to Gulp runfile

energyfile = fopen(fileenergy,'w');
for i=1:size(pot_array,1)
    
   fprintf(energyfile,'%4.0f %4.8f %4.8f\r\n',i,pot_array(i,1),pot_array(i,2));
end
fclose(energyfile);
gg
    figure
scatter3(coordjx,coordjy,coordjz,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
%view(-30,10)
axis([0 4 0 4 0 4])
set(gca,'xtick',[0:4:5])
set(gca,'ytick',[0:4:5])
set(gca,'ztick',[0:4:5])
%%
gg



    



