% % % ---------------------------------------------------------------------
% % % READ ME.
% % % This code generate the data of SYK model. Notice that this code donot
% % % plot any figure in the main text, where FIG.2 and FIG.3 are plotted in
% % % 'code3_plotting_FIG2_FIG3.m'. Instead, after running this code, 
% % % the workspace is saved as 'syk_12_data.mat' and this data is used in 
% % % code3 to plot FIG2 and FIG3. Here we use a small system size(Nf=6) to
% % % test the code, the data saved in another file is with larger size(Nf=12)
% % % ---------------------------------------------------------------------
q=4;J=1;nmax=20;
Nf=6;Ndis=10;

N=2*Nf;MajoranaS=Majorana_gen(Nf);


Dnn=zeros(nmax+1,4);
D=zeros(nmax+1,nmax+1,2);
Osize=zeros(nmax+1,2);

for idis=1:Ndis
    
%     tic
    H=H_gen(N,q,J,MajoranaS);H=full(H);
    Majorana=zeros(2^Nf,2^Nf,2*Nf);
    for i=1:N
        Majorana(:,:,i)=full(MajoranaS{i});
    end
%     toc
    
%     ['Nf=',num2str(Nf),', idis=',num2str(idis),',imode=1']
    jumprange=1;O0range=1;
    O0=O0_gen(Majorana,O0range);
    [~,tempD]=bn_gen_full(O0,H,nmax,N,Majorana,jumprange);
    D(:,:,1)=D(:,:,1)+abs(tempD);
    Dnn(:,1)=Dnn(:,1)+diag(abs(tempD));
    
    
%     ['Nf=',num2str(Nf),', idis=',num2str(idis),',imode=2']
    jumprange=1;O0range=2;
    O0=O0_gen(Majorana,O0range);
    [~,tempD]=bn_gen_full(O0,H,nmax,N,Majorana,jumprange);
    D(:,:,2)=D(:,:,2)+abs(tempD);
    Dnn(:,2)=Dnn(:,2)+diag(abs(tempD));
    
    
%     ['Nf=',num2str(Nf),', idis=',num2str(idis),',imode=3']
    jumprange=2;O0range=1;
    O0=O0_gen(Majorana,O0range);
    [~,tempDnn,tempOsize]=bn_gen(O0,H,nmax,N,Majorana,jumprange);
    Osize(:,1)=Osize(:,1)+tempOsize;
    Dnn(:,3)=Dnn(:,3)+tempDnn;  
    
    
%     ['Nf=',num2str(Nf),', idis=',num2str(idis),',imode=4']
    jumprange=2;O0range=2;
    O0=O0_gen(Majorana,O0range);
    [~,tempDnn,tempOsize]=bn_gen(O0,H,nmax,N,Majorana,jumprange);
    Osize(:,2)=Osize(:,2)+tempOsize;
    Dnn(:,4)=Dnn(:,4)+tempDnn;  
    
    
end
D=D/Ndis;Dnn=Dnn/Ndis;Osize=Osize/Ndis;

save('syk_12_data.mat');
% % % ---------------------------------------------------------------------

% % % -------This function generate the Hamitonian of SYK model, where N is
% number of Maojorana; q and J are in the definition of SYK
% model,'Majorana' is the imput of Majorana matrixes
function y=H_gen(N,q,J,Majorana)
comb=nchoosek(1:N,q);Ncomb=size(comb,1);
y=sparse(zeros(2^(N/2),2^(N/2)));
for icomb=1:Ncomb
    
    temp=Majorana{comb(icomb,1)};
    for k=2:q
        temp=temp*Majorana{comb(icomb,k)};
    end
    y=y+randn*temp;
    
end
a=sqrt( factorial(q-1) * (J^2) / (N^(q-1)) ) * ((1i)^(q/2)) ;
y=y * a;

end

% % % -------This function generate the Matrix representation of Majorana
% fermion, by means of Jorda-Wigner transformation. Nf is the fermion
% number, where Nf=2*N;
function y=Majorana_gen(Nf)
Z=[1 0;0 -1];Splus=[0 1;0 0];Sminus=[0 0;1 0];
y=cell(2*Nf,1);
for r=1:Nf
cplus=kron(kron(kronn(Z,r-1),Splus),eye(2^(Nf-r)));
cminus=kron(kron(kronn(Z,r-1),Sminus),eye(2^(Nf-r)));
y{r}=sparse((cplus+cminus)/sqrt(2));
y{r+Nf}=sparse((cplus-cminus)/(1i*sqrt(2)));
end

end

% % % -------This function generate the initial operator which is used as
% the seed of recursion process. O0range is the string length of O0, for
% example, if imput O0=2, the output is \psi_{1}*\psi_{2}
function O0=O0_gen(Majorana,O0range)
O0=Majorana(:,:,1);
if O0range>1
    for r=2:O0range
        O0=O0*Majorana(:,:,r); 
    end
end
end

% % % -------This function generate Krylov coefficient bn and full matrix
% element Dmn of the dissipator.
function [b,D]=bn_gen_full(On2,H,nmax,N,Majorana,jumprange)
dim=size(H,1);
b=zeros(nmax,1);
O=zeros(dim,dim,nmax+1);
D=zeros(nmax+1,nmax+1);
% Osize=zeros(nmax+1,1);

On2=On2/sqrt(inner(On2,On2));O(:,:,1)=On2;
On1=H*On2-On2*H;
b(1)=sqrt(inner(On1,On1));
On1=On1/b(1);
O(:,:,2)=On1;

% % % tic;'basis gen'
for n=2:nmax

        A=H*On1-On1*H-b(n-1)*On2;
        b(n)=sqrt(inner(A,A));
        On2=On1;
        On1=A/b(n);
        O(:,:,n+1)=On1;

end
% % % toc

t=size(D,1);
for I=1:(nmax+1)
% % % tic;I
    ttemp=Lind_sep(O(:,:,I),N,Majorana,jumprange);

    for J=1:t

        D(I,J)=inner(O(:,:,J),ttemp);

    end
% % % toc
end

end

% % % -------This function generate Krylov coefficeint bn and only the
% diagonal element Dnn of the disspator
function [b,Dnn,Osize]=bn_gen(On2,H,nmax,N,Majorana,jumprange)
dim=size(H,1);
b=zeros(nmax,1);
O=zeros(dim,dim,nmax+1);
Osize=zeros(nmax+1,1);
Dnn=zeros(nmax+1,1);

On2=On2/sqrt(inner(On2,On2));O(:,:,1)=On2;
On1=H*On2-On2*H;
b(1)=sqrt(inner(On1,On1));
On1=On1/b(1);
O(:,:,2)=On1;

% % % tic;'basis gen'
for n=2:nmax

        A=H*On1-On1*H-b(n-1)*On2;
        b(n)=sqrt(inner(A,A));
        On2=On1;
        On1=A/b(n);
        O(:,:,n+1)=On1;

end
% % % toc


for I=1:(nmax+1)
% % % tic;I
Dnn(I)=inner(O(:,:,I),Lind_sep(O(:,:,I),N,Majorana,jumprange));
Osize(I)=inner(O(:,:,I),Lind_sep(O(:,:,I),N,Majorana,1));
% % % toc
end

end

% % % -------This function evaluate the disspator acting on an operator O
function y=Lind_sep(O,N,Majorana,jumprange)
dim=size(O,1);
y=zeros(dim,dim);
a=[1 1i 1i 1];Hermitianfactor=a(mod(jumprange-1,4)+1);
if jumprange==1
    
    for ikai=1:N
        jump=Majorana(:,:,ikai);temp=0.5;
        y = y -  ( ctranspose(jump)*O*jump - 0.5* (temp*O+O*temp) );
    end

else
    
    jumplst=nchoosek(1:N,jumprange);
    for ijump=1:size(jumplst,1)
        jump=Majorana(:,:,jumplst(ijump,1));
        for r=2:jumprange
            jump=jump*Majorana(:,:,jumplst(ijump,r));
        end
        jump=jump*Hermitianfactor;
        temp=2^(-jumprange);
        y = y -  ( ctranspose(jump)*O*jump - 0.5* (temp*O+O*temp) );
    end
    
end
coeff=( N^(1-jumprange) ) * factorial(jumprange) * ( 2^(jumprange) );
y=y*coeff;
end

% % % -------This function generate innerproduct of operator
function y=inner(A,B)
y=trace(ctranspose(A)*B);
end

function y=kronn(A,n)
if n==0
    y=1;
elseif n==1
    y=A;
elseif n>1
    y=A;
    for k=1:(n-1)
    y=kron(y,A);
    end
    
end
end
