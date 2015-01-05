## CPHD ##
Ref: Vo B-N, Ma W-K. Analytic Implementations of the Cardinalized Probability Hypothesis Density Filter[J]. 
Signal Processing, IEEE Transactions on. 2007, 55(7) pp2559 Section IV
### 实验设置 ###

Seting basic dimension of state and observation vector. 

		nDimX = 4;    % Dimension of State vector.
        nDimZ = 2;    % Dimension of Measurement vector.

状态向量包括速度和位移，

Surveillance region.

    xmin = -1000;  % Minimum position in x axis.
    xmax = 1000;   % Maximum position in x axis.
    ymin = -1000;  % Minimum position in y axis.
    ymax = 1000;   % Maximum position in y axis.
    vol = (xmax-xmin)*(ymax-ymin);  % Volume o surveillance region.
    u = 1/vol;     % uniform distribution of clutter.

目标的活动区域是[-1000,1000]*[-1000,1000].杂波是泊松RFS,

		Clutter density.
		lambdac = 12.5e-6;%1e-6;%10e-7;%Poisson clutter spatial distribution density.

V = 4* 10 ^6 ,意味着每一帧误检的平均个数是50.

process noise and observation noise setting.

	sigmav = 10;        % standard variance of process noise.
	sigmaw = 10;       % standard variance of measurement noise.

Survival probability and detection probability.

		probsv = 1;   % Probability of target survival.
		probdt = 0.98;   % Probability of detection.
		pspar  = 0.99;   % Probability of target survival as a GM-PHD filter parameter.
		pdpar  = 0.80;   % Probability of detection as a GM-PHD filter parameter.

Threshold and maximum number of gaussian components sets

		thtrunc = 1e-5;  % threshold to truncate Gaussian component with small weight.
		thunion = 4;     % threshold to union those Gaussian component with near means.
		nGMMax = 100;    % Maximum Number of Gaussian component.
		%kspawn = 66;     % Time step at which Spawning target appears.
		nCardMax = 100;

System model and Observation model.

    F = kron(eye(2),[1 T; 0 1]);
    Q = sigmav*sigmav*kron(eye(2),[T^4/4 T^3/2; T^3/2 T^2]);
    H = kron(eye(2),[1 0]);
    R = sigmaw*sigmaw*eye(2);

详细见论文第6页。
kron函数解释：

    kron   Kronecker tensor product.
    kron(X,Y) is the Kronecker tensor product of X and Y.
    The result is a large matrix formed by taking all possible
    products between the elements of X and those of Y. For
    example, if X is 2 by 3, then kron(X,Y) is
 
       [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
         X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
 
    If either X or Y is sparse, only nonzero elements are multiplied
    in the computation, and the result is sparse.

Gaussian Mixture Model for birth target.

		GMM_birth(1).omega = 0.03;
		GMM_birth(1).mean = [-800 0 -200 0]';
		GMM_birth(1).variance = diag([10 10 10 10]);
		
		GMM_birth(2).omega = 0.03;
		GMM_birth(2).mean = [-200 0 800 0]';
		GMM_birth(2).variance = diag([10 10 10 10]);
		
		GMM_birth(3).omega = 0.03;
		GMM_birth(3).mean = [0 0 0 0]';
		GMM_birth(3).variance = diag([10 10 10 10]);
		
		GMM_birth(4).omega = 0.03;
		GMM_birth(4).mean = [400 0 -600 0]';
		GMM_birth(4).variance = diag([10 10 10 10]);

四个新生


 Generate the observation set.

        Z = [];
        for n=1:size(Xtrue,2)
            if rand < probdt,  % for probability of detection
                Z = [Z, mvnrnd((H*(Xtrue(:,n)))', R)'];
            end
        end
        nc = poissrnd(lambdac*vol);
        ZCLobsv{k+1,m} = [Z [xmin + (xmax-xmin)*rand(1,nc);...
                             ymin + (ymax-ymin)*rand(1,nc)]];


现在开始CPHD
Step1: Prediction for birth and spawning targets.

	for n=1:length(GMM_birth)
	    m = m + 1;
	    WeightPred(m) = GMM_birth(n).omega; %#ok<*AGROW>
	    MeanPred(:,m) = GMM_birth(n).mean;
	    VarPred(:,:,m) = GMM_birth(n).variance;
	end

Step2: Prediction for existing targets.

	for n=1:length(GMM)
	    m = m + 1;
	    WeightPred(m) = ps*GMM(n).omega;
	    MeanPred(:,m) = F*GMM(n).mean;
	    VarPred(:,:,m) = F*GMM(n).variance*F'+Q;
	end

对应公式25，26，27

Step3: Prediction of Cardinality distribution.  计算公式23


		ProbPred = zeros(size(Prob));
		UIdx     = triu(repmat(1:nmax+1,[nmax+1,1]),0);
		RIdx     = triu(repmat((1:nmax+1)',[1,nmax+1]),0);
		
		for n=0:nmax   
		    MatIdx = UIdx(1:n+1,:)';  % Use matrix operation for speed.
		    JJ     = MatIdx(MatIdx>0);
		    JJ     = JJ' - 1;
		    
		    MatIdx = RIdx(1:n+1,:)';
		    KK     = MatIdx(MatIdx>0);
		    KK     = KK' - 1;
		    ProbPred(n+1) = sum(factorial(JJ)./...
		                        factorial(JJ-KK)./...
		                        factorial(KK).*...
		                        (1-ps).^(JJ-KK).*...
		                        Prob(JJ+1).*...
		                        ps.^KK.*ProbBirth(n-KK+1));
		 %For Loop with slow mathod               
		    for kk=0:n  
		        hbar = sum(factorial(kk:nmax)./factorial((kk:nmax)-kk)/factorial(kk).*...
		               (1-ps).^((kk:nmax)-kk).*Prob((kk:nmax)+1));
		        ProbPred(n+1) = ProbPred(n+1) + ps^kk*ProbBirth(n-kk+1)*hbar;
		    end
		end


详细推导结合23 和 第3页的说明。

