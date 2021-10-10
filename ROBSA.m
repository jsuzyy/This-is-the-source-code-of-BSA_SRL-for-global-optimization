function [BestCost,BestValue,Best]=ROBSA(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X)
it=0;
FES=0;
FE_best=[];
DIM_RATE=1;
% Initialise the population
for i=1:n
    X(i,:)=VarMin+(VarMax-VarMin).*rand(1,nVar); % Eq. 26
end
for i=1:nPop
    val_X(i) =  fhd(X(i,:)); 
end
historical_X = repmat(VarMin, nPop, 1) + rand(nPop, nVar) .* (repmat(VarMax-VarMin, nPop, 1));
M_X=rand;
while  FES < MaxIt
    it=it+1;
    if rand<rand
        historical_X=X;
    end
    historical_X=historical_X(randperm(nPop),:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    map=ones(nPop,nVar);     %%%
    if rand<rand
        for i=1:nPop
            u=randperm(nVar);
            map(i,u(1:ceil(DIM_RATE*rand*nVar)))=0;
        end
    else
        for i=1:nPop
            map(i,randi(nVar))=0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F=3*randn;
    for i=1:nPop
        
        Xi=X(i,:)+F.*map(i,:).*(historical_X(i,:)-X(i,:));
        
        Xi = boundConstraint_absorb(Xi, VarMin, VarMax);
        val_Xi =  fhd( Xi); 
        FES = FES+1;
        if val_Xi<val_X(i)
            val_X(i) = val_Xi;
            X(i,:) = Xi;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if rand<0.3
        varmin=min(X);
        varmax=max(X);
        for i=1:nPop
            if rand<rand
                
                lamda=1+rand*rand;
            else
                lamda=1-rand*rand;
            end
            XX(i,:)=(0.5+0.5*lamda)*(varmin+varmax)-lamda*X(i,:);
            XX(i,:) = boundConstraint_absorb(XX(i,:), VarMin, VarMax);
            val_XXi(i) =  feval(fhd,(XX(i,:))',varargin{:}); 
            FES = FES+1;
        end
        XT=[X;XX];
        XCOST=[val_X val_XXi];
        [~,ind]=sort(XCOST);
        X=XT(ind(1:nPop),:);
        val_X=XCOST(ind(1:nPop));
    end
    [~,index_Best] = sort(val_X);
    BestValue=val_X(index_Best(1));
    Best=X(index_Best(1),:);
    BestCost(it)=BestValue;
  %     disp(['Iteration: ',num2str(it),'   Fmin= ',num2str(BestValue,15)]);
end
end