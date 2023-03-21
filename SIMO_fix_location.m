clear all 

N=128;%number of antennas
fc=28*10^9;%carrier frequency
lamb = 3*10^8/fc; %wavelength
d_element = lamb/2; %antenna spacing
d_ray = 2*((N-1)*d_element)^2/lamb; %Rayleigh distance 
R=0.1; %target data rate
eps = 2^R-1;
noise = 10^(-110/10); %noise power
%Ptx = 10^(-0/10);
M = 36; %number of near field users
ct=1;
Dx =4;% number of beams to be used
K=1;

%the locations of the elements
temp1 = [0:d_element:(N-1)*d_element]';
loc_array = [zeros(N,1) temp1-(N-1)*d_element/2];
snrdb = [0: 10: 30];

 
for icD =  1: length(snrdb)
    Ptx = 10^((snrdb(icD)-30)/10);
    for ict = 1 : ct
        %the location of the near-field users [5 70]
        NF_loc=[];
        sqrt_M = sqrt(M);
        step1 =  10/sqrt_M;
        temp1 = [0: step1: (sqrt_M-1)*step1];
        [X,Y] = meshgrid(temp1+5,temp1-temp1(end)/2);
        NF_loc = [X(:) Y(:)];
%         while size(NF_loc,1)<=M
%             x_loc = [d_ray*rand(1,1) sign(randn)*d_ray*rand(1,1)];
%             if sqrt(x_loc*x_loc')<d_ray & sqrt(x_loc*x_loc')>5
%                 NF_loc = [NF_loc; x_loc];
%             end
%         end

        %the channel vectors of the near-field users
        Hm = [];
        for m = 1 : M
            hmx = exp(-complex(0,1)*2*pi/lamb* sqrt(sum((NF_loc(m,:)-loc_array).^2,2)));
            Hm = [Hm hmx*3*10^8/4/pi/fc/sqrt(NF_loc(m,:)*NF_loc(m,:)')];
        end

        %the location of the far-field users [80 100]
%         K = 1; %number of near field users
%         FF_loc=[];
%         while size(FF_loc,1)<=K
%             x_loc = [70*rand(1,1) sign(randn)*70*rand(1,1)];
%             if sqrt(x_loc*x_loc')<d_ray+10 & sqrt(x_loc*x_loc')>d_ray
%                 FF_loc = [FF_loc; x_loc];
%             end
%         end
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
            costhetam =  (FF_loc(k,2))/sqrt(FF_loc(k,:)*FF_loc(k,:)');%x/r
            gkx1 = exp(-complex(0,1)*2*pi/lamb* d_element*[0:1:N-1]*costhetam );
            Gk = [Gk gkx1.'*3*10^8/4/pi/fc/sqrt(FF_loc(k,:)*FF_loc(k,:)')*...
                exp(-complex(0,1)*2*pi/lamb * sqrt(sum((FF_loc(k,:)-loc_array(1,:)).^2,2)))];
        end

        %precoding matrix based on the near-field users
        Ptemp = inv(Hm'*Hm);
        D = diag(sqrt(1./diag(Ptemp)));
        P = Hm*inv(Hm'*Hm)*D;

        %variables for the system model
        for m = 1 : M
            hm(m,1) = abs(Hm(:,m)'*P(:,m))^2; %hm
            gm(m,1) = abs(Gk'*P(:,m))^2;%gm
            Pmstar(m,1) = noise*eps/hm(m); %Pm^*
            etam(m,1) = noise + Pmstar(m)*hm(m); %etam
        end
        eta0 = noise + min(Ptx,Pmstar)'*gm;%need to cap the transmit power

        %prepare the optimization excluding infeasible beams
        ind_feasible = find(Pmstar<Ptx);
        hmz=hm(ind_feasible);gmz=gm(ind_feasible);
        Pmstarz=Pmstar(ind_feasible);etamz=etam(ind_feasible);

        %use a subset of the beams
        x_dim = min(Dx,length(Pmstarz));
        %[a, ind] = sort(gmz,'descend');%sort(min([gmz(:,k)/max(gmz(:,k)) hmz/max(hmz)],[],2),'descend');
        [a, ind] = sort(min([gmz(:,k)/max(gmz(:,k)) hmz/max(hmz)],[],2),'descend');
        S_good = ind(1:x_dim);
        hmx = hmz(S_good);
        gmx = gmz(S_good);
        Pmstarx = Pmstarz(S_good);
        etamx = etamz(S_good);

        A = []; % No other constraints
        b = [];
        Aeq = [];%[eye(M,M);-eye(M,M)];%-eye(2*Dsize);
        beq = [];%[C;zeros(M,1)];%zeros(2*Dsize,1);
        lb = [];
        ub = [];
        x0 = zeros(x_dim+1,1); %[y   xm]
        options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
        x = fmincon(@(x) -x(1),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,gmx,etamx,x_dim,eta0,hmx,Pmstarx,Ptx),options);

        [c,] = mycons(x,gmx,etamx,x_dim,eta0,hmx,Pmstarx,Ptx);
        if max(c)<=0.0000001
            rate(ict) = log2(1+x(1));
        else
            rate(ict) = 0;
        end

        if imag(rate(ict))>1
            dfd=0;
        end
    end
    rate_dx(icD) = mean(rate);
end

plot(snrdb,rate_dx)
 
function [c,ceq] = mycons(x,gmx,etamx,x_dim,eta0,hmx,Pmstarx,Ptx)
%first Dsize element of x is x, the second half is alpha^P
c(1,1) = x(1)- (sum(sqrt(gmx).*sqrt(x(2:end))))^2/eta0;
c(2:x_dim+1,1) = x(1) - x(2:end).*hmx./etamx;
c(x_dim+2:x_dim+x_dim+1,1) = x(2:end) - Ptx +  Pmstarx;
c(2*x_dim+2:3*x_dim+2,1) = -x; 
ceq = [];
end