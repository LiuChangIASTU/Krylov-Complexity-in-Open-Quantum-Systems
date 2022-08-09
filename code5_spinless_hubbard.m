function code5_spinless_hubbard
% % % ---------------------------------------------------------------------
% % % READ ME.
% % % This code generate the data of spinless fermion model. Notice that this code donot
% % % plot any figure in the main text, where FIG.2 and FIG.3 are plotted in
% % % 'code3_plotting_FIG2_FIG3.m'. Instead, after running this code, 
% % % the workspace is saved as 'fermion_13_data.mat' and this data is used in 
% % % code3 to plot FIG2 and FIG3.
% % % ---------------------------------------------------------------------
    for mode = 1:4
        spinless_fermion_krylov(mode);    %use system size N = 8 to test
    end
end

function spinless_fermion_krylov(mode)
% 1d spinless model with NN and NNN hopping, NN and NNN interaction
%mode is parameter for initial krylov basis and jump operator M
%mode = 1 O = n_1               M_i = 2n_i
%mode = 2 O = n_1               M_i = f_i+f_i^\dagger
%mode = 3 O = f_1+f_1^\dagger   M_i = 2n_i
%mode = 4 O = f_1+f_1^\dagger   M_i = f_i+f_i^\dagger

    N = 8;     %system size
    t1 = 1;     % NN hopping
    t2 = 0.2;   % NNN hopping
    u1 = 0.6;   % NN interaction
    u2 = 0.1;   % NNN interaction
    gamma = 1;  %dissipation strength
    dim = 50;   %krylov subspace dimension
    
    %construct fermion operator
    psif = cell(N,1);       %complex fermion annihilation operator
    c = sparse([0,1;0,0]);  %single fermion annihilation operator
    cd = c';                %creation operator
    sz = sparse([1,0;0,-1]);%pauli sigma z
    for count = 1:N
        if count ==1
            jw = 1;
        else
            jw = kron(sz,jw);% jordan wigner factor
        end
        psif{count,1} = kron(kron(jw,c),eye(2^(N-count)));
    end
    
    %construct Hmailtonian
    H = sparse(zeros(2^N,2^N));
    for count = 1:N
        f1 = psif{count,1};
        f2 = psif{mod(count,N)+1,1};
        f3 = psif{mod(count+1,N)+1,1};
        
        H = H-t1*(f1'*f2+f2'*f1);
        H = H-t2*(f1'*f3+f3'*f1);
        
        H = H+u1*(f1'*f1*f2'*f2);
        H = H+u2*(f1'*f1*f3'*f3);
    end 
    
    %construct kryloc subspace
    if mode == 1
        O_ini =psif{1,1}'*psif{1,1};
    elseif mode == 2
        O_ini =psif{1,1}'*psif{1,1};
    elseif mode == 3
        O_ini =psif{1,1}+psif{1,1}';
    elseif mode == 4
        O_ini =psif{1,1}+psif{1,1}';
    else
        disp('please choose mode 1-4');
        return
    end
    %krylov basis: k_basis
    %Lanczos coffecient: bn
    %dimension = dim
    [bn,k_basis] = krylov_basis(H,O_ini,dim);
    
    
    %calculate dmn, dn
    Dmn = zeros(dim,dim); %D matrix

    % if calculate operator size of krylov basis,construct these variables
    %otoc = zeros(dim,1);  %otoc
    %opsize_k = zeros(dim,1);  %operator size
    
    for count1 = 1:dim 
        Dtemp = sparse(zeros(2^N,2^N));
        for num = 1:N
            f_temp = psif{num,1};
            if mode == 1 || mode == 3
                M_temp = 2*f_temp'*f_temp;  %jump operator is 2n_i   boson
            elseif mode == 2 || mode == 4
                M_temp = f_temp+f_temp'; %jump operator is f+f^\dagger   fermion
            end        
            kb_temp = full(k_basis{count1,1});%krylov basis
             
            if mode == 4   %fermion O fermion M
                Dtemp = Dtemp+0.5i*gamma*M_temp'*M_temp*kb_temp+0.5i*gamma*kb_temp*M_temp'*M_temp+1i*gamma*M_temp'*kb_temp*M_temp;
            else
                Dtemp = Dtemp+0.5i*gamma*M_temp'*M_temp*kb_temp+0.5i*gamma*kb_temp*M_temp'*M_temp-1i*gamma*M_temp'*kb_temp*M_temp;
            end
        end     


        for count2 = 1:dim
            kb_temp = full(k_basis{count2,1});
            Dmn(count1,count2) = trace(kb_temp'*Dtemp);
        end
        
        
        %operator size
        %{
        for count_s = 1:N
            kb_temp = full(k_basis{count1,1});
            psim1 = psif{count_s,1}+psif{count_s,1}';     %maorana fermion
            psim2 = (psif{count_s,1}-psif{count_s,1}')/(1i);
                        
            otoc(count1) = otoc(count1)+trace(kb_temp'*psim1*kb_temp*psim1);
            otoc(count1) = otoc(count1)+trace(kb_temp'*psim2*kb_temp*psim2);
        end
        if mode == 1 || mode == 2
            opsize_k(count1) = (-real(otoc(count1))+2*N)/2;  %boson
        else
            opsize_k(count1) = (real(otoc(count1))+2*N)/2;  %fermion
        end
        %}

    end
    %absolute value
    Dmn = abs(Dmn);
    dn = diag(Dmn);

    
    % save N=13 data 
    % first need choose mode = 1 to create fermion_13_data.mat
    if mode ==1
        bb = dn;
        save('fermion_13_data.mat','bb');
    elseif mode==2
        bf = dn;
        save('fermion_13_data.mat','bf','-append');
        save('fermion_13_data.mat','Dmn','-append');
    elseif mode==3
        fb = dn;
        save('fermion_13_data.mat','fb','-append');
    elseif mode==4
        ff = dn;
        save('fermion_13_data.mat','ff','-append');
    end
    
    
end





function [bn,k_basis] = krylov_basis(H,O_ini,D) %O_ini is O_0, first krylov bsis %basis number cut by D
    d = length(H);      %hilbert space dim
    bn = zeros(D,1);
    %bn(1)=b0 =0
    k_basis = cell(D,1);
    %construct k_basis and bn
    k_basis{1,1} = O_ini/sqrt(trace(O_ini'*O_ini)); %normalization
    k_basis{2,1} = H*k_basis{1,1}-k_basis{1,1}*H;   %A_1
    bn(2) = sqrt(trace(k_basis{2,1}'*k_basis{2,1}));%b1
    k_basis{2,1} = k_basis{2,1}/bn(2);              %normalization
    for count = 3:D
        
        kb1 = full(k_basis{count-1,1});
        kb2 = full(k_basis{count-2,1});
        
        temp = H*kb1-kb1*H-bn(count-1)*kb2;
        bn(count) = sqrt(trace(temp'*temp));
        k_basis{count,1} = sparse(temp/bn(count));
    end

end