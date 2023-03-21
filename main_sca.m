clear all 
figure
N=64;%number of antennas
fc=28*10^9;%carrier frequency
lamb = 3*10^8/fc; %wavelength
d_element = lamb/2; %antenna spacing
d_ray = 2*((N-1)*d_element)^2/lamb; %Rayleigh distance 
d_ray_128 = 86.4054;%this is the Rayleigh distance for the case of N=128
d_ray_64 = 21.2625;
R=0.1; %target data rate
eps = 2^R-1;
noise = 10^(-110/10); %noise power
Ptx = 10^(-0/10);
M = 36; %number of near field users
K=1;%number of far field users
ct=[100 400 1000 1000];
Dx =1;% number of beams to be used
eps_BB = 0.01;
max_itr=60;

%the locations of the elements
temp1 = [0:d_element:(N-1)*d_element]';
loc_array = [zeros(N,1) temp1-(N-1)*d_element/2];

snrdb = [0: 10: 30];


for icD =   1:  length(snrdb) 
    Ptx = 10^((snrdb(icD)-30)/10);
    num_event(icD) = 0; 
    for ict = 1 : ct(icD)
        %the location of the near-field users [5 70]
        NF_loc=[];
        while size(NF_loc,1)<=M
            x_loc = [d_ray_64*rand(1,1) sign(randn)*d_ray_64*rand(1,1)];
            if sqrt(x_loc*x_loc')<d_ray_64 & sqrt(x_loc*x_loc')>5
                NF_loc = [NF_loc; x_loc];
            end
        end

        %the channel vectors of the near-field users
        Hm = [];
        for m = 1 : M
            hmx1 = exp(-complex(0,1)*2*pi/lamb* sqrt(sum((NF_loc(m,:)-loc_array).^2,2)));
            Hm = [Hm hmx1*3*10^8/4/pi/fc/sqrt(NF_loc(m,:)*NF_loc(m,:)')];
        end

        %the location of the far-field users [80 100]
        FF_loc=[];
        while size(FF_loc,1)<=K
            x_loc = [70*rand(1,1) sign(randn)*70*rand(1,1)];
            if sqrt(x_loc*x_loc')<d_ray_128+10 & sqrt(x_loc*x_loc')>d_ray_128
                FF_loc = [FF_loc; x_loc];
            end
        end

        %channel vectors of the far-field users
        Gk = [];
        ggk = [];
        for k = 1 : K
            costhetam =  (FF_loc(k,2))/sqrt(FF_loc(k,:)*FF_loc(k,:)');%y/r
            gkx1 = exp(-complex(0,1)*2*pi/lamb* d_element*[0:1:N-1]*costhetam );
            Gk = [Gk gkx1.'*3*10^8/4/pi/fc/sqrt(FF_loc(k,:)*FF_loc(k,:)')*...
                exp(-complex(0,1)*2*pi/lamb * sqrt(sum((FF_loc(k,:)-loc_array(1,:)).^2,2)))];
        end

        %precoding matrix based on the near-field users
        Ptemp = inv(Hm'*Hm);
        D = diag(sqrt(1./diag(Ptemp)));
        P = Hm*inv(Hm'*Hm)*D;
        Gk_P = P'*Gk; %tilde{g}_k


        %variables for the system model
        for m = 1 : M
            hm(m,1) = abs(Hm(:,m)'*P(:,m))^2; %hm
            gm(m,:) = abs(Gk'*P(:,m)).^2.';%K column, M row, |gk'*pm|^2
            Pmstar(m,1) = noise*eps/hm(m); %Pm^*
            etam(m,1) = (noise + Pmstar(m)*hm(m))/hm(m); %etam
        end
        eta0k = noise +  min(Ptx, Pmstar)'*gm;%need to cap the transmit power
        gktilde = P'*Gk;%M rows, K column, no abs

        %prepare the optimization excluding infeasible beams
        ind_feasible = find(Pmstar<Ptx);
        hmz=hm(ind_feasible);gmz=gm(ind_feasible,:); 
        Pmstarz=Pmstar(ind_feasible);etamz=etam(ind_feasible);
        gktildez = gktilde(ind_feasible,:);

        %scheduling
        num_beam = length(Pmstarz); %number of qualified beams
        if num_beam<=(K-1)*Dx
            num_event(icD) = num_event(icD)+1
            continue;
        end
        hmz_temp= hmz; gmz_temp= gmz; Pmstarz_temp= Pmstarz; 
        etamz_temp= etamz; gktildez_temp= gktildez; %this is to back up gmz
        for k =1 : K %we have to use arrays, since some users might get less beams than others
            x_dim(k) = min(Dx,size(gmz_temp,1)); %the number of beams to be used by user k
            %[a, ind] = sort(gmz_temp(:,k),'descend');
            [a, ind] = sort(min([gmz_temp(:,k)/max(gmz_temp(:,k)) hmz_temp/max(hmz_temp)],[],2),'descend');
            S_good{k} = ind(1:x_dim(k)); %indeces of the beams to by used by the user
            hmx{k} = hmz_temp(S_good{k});
            gmx{k} = gmz_temp(S_good{k},:);
            Pmstarx{k} = Pmstarz_temp(S_good{k});
            etamx{k} = etamz_temp(S_good{k});
            gktildezx{k} = gktildez_temp(S_good{k});
            gmz_temp(S_good{k},:) = [];%by setting them zero, avoid the beams to be used by other users
            hmz_temp(S_good{k},:) = []; 
            Pmstarz_temp(S_good{k},:) = []; 
            etamz_temp(S_good{k},:) = [];
            gktildez_temp(S_good{k},:) = [];
        end

        %initilization
        x_taylor=[];
        for k = 1 : K
            fk{k} = max(0,(Ptx-Pmstarx{k})).*ones(x_dim(k),1);
        end
        for k = 1 : K
            sum1=0;
            for i = 1 : K
                if i==k
                    continue;
                else
                    sum1 = sum1 + abs(gktildez(S_good{i},k)'*fk{i})^2; 
                    %the k-th column for user k, S_good{i} are the interferring beams by user i
                end
            end
            temp1 = [abs(gktildezx{k}'*fk{k})^2/(eta0k(k)+ sum1)];
            mum = (noise+sum(Pmstarx{k}.*hmx{k}))./hmx{k};
            x_taylor(k,1) = min([temp1 ;abs(fk{k}).^2./mum]);
        end
        for k = 1 : K
            x_taylor = [x_taylor;real(fk{k});imag(fk{k})];
        end
        rate_ini = sum(log2(1+x_taylor(1:K)));

        A = []; % No other constraints
        b = [];
        Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
        beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
        lb = [];
        ub = [];
        x_taylor_backup = x_taylor;

        iter = 1; error = 100000;
        rate_feasible=0;
        while iter< max_itr & error>0.001
            temp2=temp1;
            x0 =  x_taylor;%zeros(length(x_taylor),1); %[x_1 ... x_K real(f_1) image(f_1)... ... f_K]
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            x = fmincon(@(x) -sum(log2(1+x(1:K))),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good),options);
            %[cx,] =mycons(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good);
            
            temp1 = sum(log2(1+x(1:K)));
            
            %save the feasible solutions
            [c2,]=feasibility(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good); %y is matrix 
            if max(c2)<0.00000000001
                rate_feasible = max([temp1 rate_feasible rate_ini]);
            end
            %[max(cx) temp1-temp2]
            zzz=x_taylor;
            error = abs(temp1-temp2)/abs(temp1);%sum(abs(x-x0));
            x_taylor = x;
            iter = iter+1;
 
%             if (error<0.01&iter>20)
%                 break;
%             end
        end

        rate(ict) = rate_feasible;
        
        rate_r(ict) = rate_ini; 
        % test%%%%%%%%%%%%
        %[c,c3]=feasibility(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx) %y is matrix
%         [c2,]=feasibility(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good); %y is matrix 
%         if max(c2)>0.00000000001
%             rate(ict) = 0;
%         else
%             rate(ict) = sum(log2(1+x(1:K)));
%             rateall(ict,:) = log2(1+x(1:K));
%         end

    end
    rate_dx(icD) = mean(rate)
    rate_ran(icD) = mean(rate_r)
end

plot(snrdb,rate_ran,snrdb,rate_dx)
 

 
 
function [c,ceq] = mycons(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) %y is matrix
%first Dsize element of x is x, the second half is alpha^P
%extra fk
for k = 1 : K
    fk{k} = x(K+2*sum(x_dim(1:k-1))+1:K+2*sum(x_dim(1:k-1))+x_dim(k))...
        +1i*x(K+2*sum(x_dim(1:k-1))+x_dim(k)+1:K+2*sum(x_dim(1:k)));
    fk0{k} = x_taylor(K+2*sum(x_dim(1:k-1))+1:K+2*sum(x_dim(1:k-1))+x_dim(k))...
        +1i*x_taylor(K+2*sum(x_dim(1:k-1))+x_dim(k)+1:K+2*sum(x_dim(1:k)));
end
for k = 1 : K
    sum1=0;
    for i = 1 : K
        if i == k
            continue;
        else
            sum1 =  sum1 + abs(gktildez(S_good{i},k)'*fk{i})^2;
        end
    end
    %break it to real and imag parts
    fbar = [real(fk{k});imag(fk{k})]; ftilde = [fbar; x(k)];
    f0bar = [real(fk0{k});imag(fk0{k})];f0tilde = [f0bar; x_taylor(k)];
    ghat = [real(gktildezx{k});imag(gktildezx{k})];
    gcheck = [ -imag(gktildezx{k});real(gktildezx{k})];
    
    
    diffe = [2* (ghat*ghat'+gcheck*gcheck') * f0bar./x_taylor(k); -abs(gktildezx{k}'*fk0{k})^2/x_taylor(k)^2 ];
    c(k,1) = 1 + sum1/eta0k(k) - abs(gktildezx{k}'*fk0{k})^2/x_taylor(k)/eta0k(k)...
         -2*real( diffe'*(ftilde - f0tilde ) )/eta0k(k);   
%    c(k,1) = fk{k}'*fk{k}-100000000000;
end 
for k = 1 : K
    temp1 = - Ptx +Pmstarx{k};
    c(K+1+(k-1)*x_dim(k):K+k*x_dim(k),1) =  abs(fk{k}).^2 +temp1 ;%- Ptx +Pmstarx{k}];
end
for k = 1 : K
    mum = (noise+sum(Pmstarx{k}.*hmx{k}))./hmx{k};
    c = [c; x(k)-abs(fk0{k}).^2./mum-4*real( conj(fk0{k}) .*(fk{k}-fk0{k}))./mum];
end
c = [c;-x(1:K)];

ceq = [];
end

function [c,ceq] = feasibility(x,x_taylor,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) %y is matrix
%first Dsize element of x is x, the second half is alpha^P
%extra fk
for k = 1 : K
    fk{k} = x(K+2*sum(x_dim(1:k-1))+1:K+2*sum(x_dim(1:k-1))+x_dim(k))...
        +1i*x(K+2*sum(x_dim(1:k-1))+x_dim(k)+1:K+2*sum(x_dim(1:k)));
end
for k = 1 : K
    sum1=0;
    for i = 1 : K
        if i == k
            continue;
        else
            sum1 =  sum1 + abs(gktildez(S_good{i},k)'*fk{i})^2;  
        end
    end
    c(k,1) = eta0k(k) + sum1 - abs(gktildezx{k}'*fk{k})^2/x(k) ;   
%    c(k,1) = fk{k}'*fk{k}-100000000000;
end 
for k = 1 : K
    temp1 = - Ptx +Pmstarx{k};
    c(K+1+(k-1)*x_dim(k):K+k*x_dim(k),1) =  abs(fk{k}).^2 +temp1 ;%- Ptx +Pmstarx{k}];
end
for k = 1 : K
    mum = (noise+sum(Pmstarx{k}.*hmx{k}))./hmx{k};
    c = [c; x(k)*mum-abs(fk{k}).^2];
end


ceq = [];
end  
 