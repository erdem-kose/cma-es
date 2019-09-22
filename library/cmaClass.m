classdef cmaClass
    % developed on PureCMA-ES
    % http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#matlab
    
    properties
        %input parameters
        ChrLength=3; % Chromosome length
        Lambda=100; % population size, offspring number
        Mu=100; % number of parents/points for recombination
        FunEvl=1000; % maximum function evaluation count
        Sigma=0.1; %search variance
        UpLim=1; %Chromosome top limits, (1,ChrLength) vector or single value must be uset
        LwLim=0; %Chromosome bottom limits, (1,ChrLength) vector or single value must be uset
         
        %output indicators
        solBest;
        solCand;
        CntGen=0; % final generation count
        CntEvl=0; % final evaluation count
    end
    
    methods
        function cmaObj = optimize(cmaObj, costFunc)
            %% --------------------  Initialization -------------------------------- 
            
            % Random population initialization
            Population=(cmaObj.UpLim-cmaObj.LwLim).*cmaObj.chaoRnd(cmaObj.ChrLength,cmaObj.Lambda)+cmaObj.LwLim;

            % Initial solution point
            costVals = costFunc(Population);
            [~, midx] = min(costVals);
            xmean = Population(:,midx);  

            % Strategy parameter setting: Selection  
            mu = cmaObj.Mu;
            weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
            mu = floor(mu);        
            weights = weights/sum(weights);     % normalize recombination weights array
            mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

            % Strategy parameter setting: Adaptation
            cc = (4 + mueff/cmaObj.ChrLength) / (cmaObj.ChrLength+4 + 2*mueff/cmaObj.ChrLength); % time constant for cumulation for C
            cs = (mueff+2) / (cmaObj.ChrLength+mueff+5);  % t-const for cumulation for cma_Pars.Sigma control
            c1 = 2 / ((cmaObj.ChrLength+1.3)^2+mueff);    % learning rate for rank-one update of C
            cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((cmaObj.ChrLength+2)^2+mueff));  % and for rank-mu update
            damps = 1 + 2*max(0, sqrt((mueff-1)/(cmaObj.ChrLength+1))-1) + cs; % damping for cma_Pars.Sigma 
                                                              % usually close to 1
            % Initialize dynamic (internal) strategy parameters and constants
            pc = zeros(cmaObj.ChrLength,1); ps = zeros(cmaObj.ChrLength,1);   % evolution paths for C and cma_Pars.Sigma
            B = eye(cmaObj.ChrLength,cmaObj.ChrLength);                       % B defines the coordinate system
            D = ones(cmaObj.ChrLength,1);                      % diagonal D defines the scaling
            C = B * diag(D.^2) * B';            % covariance matrix C
            invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
            eigeneval = 0;                      % track update of B and D
            chiN=cmaObj.ChrLength^0.5*(1-1/(4*cmaObj.ChrLength)+1/(21*cmaObj.ChrLength^2));  % expectation of 
                                              %   ||cmaObj.ChrLength(0,I)|| == norm(randn(cmaObj.ChrLength,1))
            % -------------------- Generation Loop --------------------------------
            cmaObj.CntEvl = 0;  % the next 40 lines contain the 20 lines of interesting code 
            cmaObj.CntGen = 0;
            
            cmaObj.solCand.Chromosome=[];
            cmaObj.solCand.Cost=0;
            bestSolTemp(1,1:uint64(floor(cmaObj.FunEvl/cmaObj.Lambda)))=cmaObj.solCand;
            cmaObj.solCand=bestSolTemp;
            clear bestSolTemp;
            
            rng('shuffle','combRecursive');
            while cmaObj.CntEvl < cmaObj.FunEvl
                cmaObj.CntGen = cmaObj.CntGen + 1;
                % Generate and evaluate cmaObj.Lambda offspring
                arx = zeros(cmaObj.Lambda,cmaObj.ChrLength);
                for k=1:cmaObj.Lambda
                    newoff = xmean + cmaObj.Sigma .* (B * (D .* randn(cmaObj.ChrLength,1))); % m + sig * Normal(0,C)
                    newoff(newoff>cmaObj.UpLim)=cmaObj.UpLim(newoff>cmaObj.UpLim);
                    newoff(newoff<cmaObj.LwLim)=cmaObj.LwLim(newoff<cmaObj.LwLim);
                    
                    if max(max(isnan(newoff)))==1
                        return;
                    elseif max(max(isinf(newoff)))==1
                        return;
                    end
                    arx(k,:)=newoff;
                    cmaObj.CntEvl = cmaObj.CntEvl+1;
                end
                arcost = costFunc(arx');

                % Sort by fitness and compute weighted mean into xmean
                [arcost, arindex] = sort(arcost);  % minimization
                xold = xmean;
                xmean = arx(arindex(1:mu),:)' * weights;  % recombination, new mean value
                cmaObj.solCand(cmaObj.CntGen).Chromosome = arx(arindex(1),:)';
                cmaObj.solCand(cmaObj.CntGen).Cost = arcost(1);

                % Cumulation: Update evolution paths
                ps = (1-cs) * ps ... 
                      + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / cmaObj.Sigma; 
                hsig = sum(ps.^2)/(1-(1-cs)^(2*cmaObj.CntEvl/cmaObj.Lambda))/cmaObj.ChrLength < 2 + 4/(cmaObj.ChrLength+1);
                pc = (1-cc) * pc ...
                      + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / cmaObj.Sigma; 

                % Adapt covariance matrix C
                artmp = (1/cmaObj.Sigma) * (arx(arindex(1:mu),:)' - repmat(xold,1,mu));  % mu difference vectors
                C = (1-c1-cmu) * C ...                   % regard old matrix  
                     + c1 * (pc * pc' ...                % plus rank one update
                             + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
                     + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 

                % Adapt step size cma_Pars.Sigma
                cmaObj.Sigma = cmaObj.Sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 

                % Update B and D from C
                if cmaObj.CntEvl - eigeneval > cmaObj.Lambda/(c1+cmu)/cmaObj.ChrLength/10  % to achieve O(cmaObj.ChrLength^2)
                  eigeneval = cmaObj.CntEvl;
                  C = triu(C) + triu(C,1)'; % enforce symmetry
                  if max(max(isnan(C)))==1
                      return;
                  elseif max(max(isinf(C)))==1
                      return;
                  end
                  [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
                  D = sqrt(diag(D));        % D contains standard deviations now
                  invsqrtC = B * diag(D.^-1) * B';
                end
            end
            solCandMat=squeeze(cell2mat(struct2cell(cmaObj.solCand)));
            [~,solBestInd]=min(solCandMat(4,:));
            cmaObj.solBest=solCandMat(1:cmaObj.ChrLength,solBestInd);
        end % end function
        function x=chaoRnd(~,m,n)
            %https://en.wikipedia.org/wiki/Logistic_map
            
            %Pomeau–Manneville scenario: chaotic in (3.56995,3.82843)
            PMs_up_bnd=3.82843; PMs_lw_bnd=3.56995;
            
            rng('shuffle','combRecursive');
            x=zeros(m,n);
            x(1,1)=rand();
            
            r=(PMs_up_bnd-PMs_lw_bnd).*rand(m,1)+PMs_lw_bnd;

            for i=2:m
                x(i,1)=r(i).*x(i-1,1).*(1-x(i-1,1));
            end
            for i=2:n
                x(:,i)=r.*x(:,i-1).*(1-x(:,i-1));
            end
        end
    end
end

