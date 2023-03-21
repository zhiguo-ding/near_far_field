clear all 

N=64;%number of antennas
fc=28*10^9;%carrier frequency
lamb = 3*10^8/fc; %wavelength
d_element = lamb/2; %antenna spacing
d_ray = 2*((N-1)*d_element)^2/lamb; %Rayleigh distance 
R=0.1; %target data rate
eps = 2^R-1;
noise = 10^(-110/10); %noise power
Ptx = 1;
M = 36; %number of near field users
K=4;%number of far field users
ct=1;
Dx =1;% number of beams to be used
eps_BB = 0.001;
max_itr=1000; % for K=4, 1000 is needed

%the locations of the elements
temp1 = [0:d_element:(N-1)*d_element]';
loc_array = [zeros(N,1) temp1-(N-1)*d_element/2];
snrdb = [0: 10: 30];

 
for icD =   length(snrdb) 
    Ptx = 10^((snrdb(icD)-30)/10);
    for ict = 1 : ct
        NF_loc=[];
        sqrt_M = sqrt(M);
        step1 = 10/sqrt_M;
        temp1 = [0: step1: (sqrt_M-1)*step1];
        [X,Y] = meshgrid(temp1+5,temp1-temp1(end)/2);
        NF_loc = [X(:) Y(:)];

        %the channel vectors of the near-field users
        Hm = [];
        for m = 1 : M
            hmx1 = exp(-complex(0,1)*2*pi/lamb* sqrt(sum((NF_loc(m,:)-loc_array).^2,2)));
            Hm = [Hm hmx1*3*10^8/4/pi/fc/sqrt(NF_loc(m,:)*NF_loc(m,:)')];
        end

        %the location of the far-field users [80 100]
        FF_loc=[];
        step_theta = pi/(K+1);
        ff_theta = [1:1:K]*step_theta-pi/2;%[-pi/2 :step_theta: pi/2];ff_theta=ff_theta(1:K);
        d = 90;%d_ray+10;
        for k = 1 : K
            FF_loc = [FF_loc; [d*cos(ff_theta(k)) d*sin(ff_theta(k))]];
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

        hm=[];gm=[];Pmstar=[];etam=[];
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
            mum = (noise+sum(Pmstarx{k}.*hmx{k}))./hmx{k};
            x_ini(k,1) = min(abs(gktildezx{k})^2*fk{k}/eta0k(k) , fk{k}/mum );
        end

        Bset=[zeros(K,1) x_ini];
        if feasibility(Bset(:,1), x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) 
            Upk = -sum(log2(1+Bset(:,1)));
            Lowk = -sum(log2(1+Bset(:,2)));
        else
            Upk = 0;Lowk = 0;
        end
        Bset = [Bset; Lowk Upk];%odd position is for lower bound
        % the last row of Bset is for low and upper bound values

        bb_k=1;
        low_ind = 1; %index of the set whose lower bound is highest
        while (Upk(end)-Lowk(end)>eps_BB) & (bb_k<max_itr)
            vec_low = Bset(end,1:2:end);% lower bound of each set, odd positions
            [temp,low_ind] = min(vec_low);            
            B_temp = Bset(1:end-1,2*(low_ind-1)+1:2*(low_ind-1)+2);%find the corresponding set 
            Bsetx = Bset(:,2*(low_ind-1)+1:2*(low_ind-1)+2);
            Bset(:,2*(low_ind-1)+1:2*(low_ind-1)+2) = [];%remove the set
            length_temp = B_temp(1:end,2)-B_temp(1:end,1);%the last row is not used
            [tempbb, ind_edge] = max(length_temp);%find the longest edge

            Bset_low_half = B_temp;  
            Bset_low_half(ind_edge, 2) = sum(Bset_low_half(ind_edge, :))/2;
            if feasibility(Bset_low_half(:,1), x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) 
                Lowk1 = -sum(log2(1+Bset_low_half(:,2)));
                Upk1 = -sum(log2(1+Bset_low_half(:,1)));
            else
                Lowk1 = 0;
                Upk1 = 0;
            end
            Bset_up_half = B_temp;  
            Bset_up_half(ind_edge, 1) = sum(Bset_up_half(ind_edge, :))/2;
            if feasibility(Bset_up_half(:,1), x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) 
                Lowk2 = -sum(log2(1+Bset_up_half(:,2)));
                Upk2 = -sum(log2(1+Bset_up_half(:,1)));
            else
                Lowk2 = 0;
                Upk2 = 0;
            end        

            Bset = [Bset [Bset_low_half; Lowk1 Upk1]]; %add this new set
            Bset = [Bset [Bset_up_half; Lowk2 Upk2]]; %add this new set
            %Bounding 
            Lowknew = min(Bset(end,1:2:end));%the new lower bound
            Lowk = [Lowk Lowknew]; %  need to remember this index
            Upnew = min(Bset(end,2:2:end));%the new upper bound
            Upk = [Upk Upnew]; % need to remember this index
            %[Lowk; Upk]
            %prunning
            vec_pru = Bset(end,1:2:end);% lower bound of each set, odd positions
            ind_pru = find(vec_pru>Upk(end));
            if ~isempty(ind_pru)
                Bset(:,[2*(ind_pru-1)+1 2*(ind_pru-1)+2])=[];
            end
            bb_k = bb_k+1;

            %[Upk;Lowk]
        end

        %recover the x
        [tempbb,opt_ind] = min(Bset(end,2:2:end));%find which box contributes the upper bound
        Boptimal = Bset(1:end-1,2*(opt_ind-1)+1:2*(opt_ind-1)+2);%extract the box

        x= Boptimal(:,1);
        rate(ict) = sum(log2(1+x));

    end

    rate_dx(icD) = mean(rate);
end

%plot(Dx_vec,rate_dx)
% plot(K_vec,rate_dx)
plot(snrdb,rate_dx)
 

 
%bisection methods to find the projection on G 
function [temp,y] = feasibility(x, x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) 
    A = []; % No other constraints
    b = [];
    Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
    beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
    lb = [];
    ub = [];    
    
    %x contains [f1... fK] 
    f0 = zeros(K,1);%zeros(length(x_taylor),1); %[x_1 ... x_K real(f_1) image(f_1)... ... f_K]
    options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
    f = fmincon(@(y) 1,f0,A,b,Aeq,beq,lb,ub,@(f) mycons(f,x,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good),options);
     
 
    [c,] = mycons(f,x,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) ;       
    if max(c)<=0.000001
        temp = 1; %feasible
    else
        temp = 0; %infeasible
    end
end

 

 
function [c,ceq] = mycons(f,x,x_dim,eta0k,gktildezx,Pmstarx,hmx,K,noise,Ptx,gktildez, S_good) 
 
for k = 1 : K
    sum1=0;
    for i = 1 : K
        if i == k
            continue;
        else
            sum1 =  sum1 + abs(gktildez(S_good{i},k)'*f(i))^2;
        end
    end    
    
    c(k,1) = x(k) + x(k)*sum1/eta0k(k) - abs(gktildezx{k})^2*f(k)/eta0k(k);   
end 
for k = 1 : K
    temp1 = - Ptx +Pmstarx{k};
    c(K+1+(k-1)*x_dim(k):K+k*x_dim(k),1) =  f(k) +temp1 ;%- Ptx +Pmstarx{k}];
end
for k = 1 : K
    mum = (noise+sum(Pmstarx{k}.*hmx{k}))./hmx{k};
    c = [c; x(k)-f(k)/mum];
end
c = [c;-f(1:K)];

ceq = [];
end
 
 