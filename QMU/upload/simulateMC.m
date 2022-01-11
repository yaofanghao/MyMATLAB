function [funcVal, stdev, CI, MCfuncVals, MCsamples]= simulateMC(func,dists,N,varargin)
% simulateMC() is a tool for Monte Carlo simulation. It draws many, many 
% samples from user-specified distributions and combines the samples through 
% an arbitrary function. 

% INPUTS:
%   func            Function handle of the function to combine samples through
%                   (func must be vectorizable so use .*, ./ and .^ instead of *, /, and ^)
%   dists           Vertical cell array of cell arrays containing params for each distribution
%   N               Number of samples to generate from each distribution
%   'CI',value      (Optional) Value to use for pseudo CI threshold. Default is 0.68.
%   'hist',nBins    (Optional) Plot histogram of the function outputs. nBins value is optional
%   'varHist',nBins (Optional) Plot histograms of each sample. nBins valueis optional.
%   'mean'|'median'|'max'|'min'
%                   (Optional) Specifies what value to return for 'funcVal'. 
%                    Default is the center of the CI inteval.
%
%   For each cell array inside 'dists' the syntax is:
%       {'DistributionType',param_1,...,paramN[,'truncate',lowerBound,upperBound]}
%   where the truncation arguments inside the square brackets are optional.
%
%   A list of supported distributions and their parameters can be found here:
%       https://www.mathworks.com/help/stats/supported-distributions.html
%   'Custom', 'Bootstrap', 'BootstrapMean', 'Fit' and 'Corr' are  additional 
%   options for 'DistributionType'. To see more about them examine the 
%   relevant comments in the code below and look in 'simulateMC_examples.m'.

% OUTPUTS: 
%   funcVal     Calculated function value. Defaults to center of CI. Can be mean, median, or max/min.
%   stdev       Standard deviation of the function outputs
%   CI          Vector containing the lowerbound and upperbound of the pseudo CI interval.
%               (this is the shortest interval containing some percent of the function outputs)
%   MCfuncVals  Vector of the function outputs.
%   MCsamples   Array containing the samples for each distribution.

% %EXAMPLE 1
% %Combine samples from X and Y through the function 'f'. 
% %'X' is a normal distribution with params [mean=5,sigma=1]. 
% %'Y' is a uniform distribution with params [lowerbound=6,upperbound=7].
% f=@(X,Y)5.*X+Y;
% dists = {{'Normal',5,1};
%          {'Uniform',6,7}};
% n=100000;
% [funcVal,stdev,CI] = simulateMC(f,dists,n)

% %EXAMPLE 2
% %Combine samples from W, X, Y, and Z through the function 'f'.
% %'W' is a custom data sample provided by the user.
% %'X' is a truncated normal distribution with params [mean=4,sigma=0.2,lowerB=3.5,upperB=4.5]. 
% %'Y' is a triangular distribution that is correlated with Z. Params are [lowerB=3.3,peak=3.5,upperB=3.7].
% %'Z' is a lognormal distribution that is correlated with Y. Params are [mean=1,sigma=0.05].
% %The correlation betwen Y and Z is specified by the Pearson correlation matrix 'corrMat'.
%
% %Set a pseudo CI of 0.95, plot histograms of all the data, and return the mean for funcVal.
%
% n=100000;
% W_data = normrnd(5,1,n,1);
% f=@(W,X,Y,Z)X.^2./(Y.*Z)+W;
% corrMat = [1 0.5; 0.5 1];
% dists = {{'Custom',W_data}
%          {'Normal',4,0.2,'trunc',3.5,4.5};
%          {'Corr',{{'Triangular',3.3,3.5,3.7};
%                   {'Lognormal',1,0.05}},corrMat}};
% [val,stdev,CI,MCfuncVals,MCsamples] = simulateMC(f,dists,n,'CI',0.95,'hist','varHist',65,'mean');

% % Note: Correlated or truncated samples are not supported for distribution
% % types of 'Bootstrap', 'Bootstrapmean', 'Custom', or 'Fit'.

% See 'simulateMC_examples.m' for more detailed examples.
% Joe Klebba, 3/2021.

%% Parse varargin to handle optional arguments
%Init defaults
if (~exist('N', 'var'))
   N=100000; 
end
confInterval = 0.68; 
meanFlag=0; medianFlag=0; maxFlag=0; minFlag=0;
histFlag=0; histBins=40;
varHistFlag = 0; varHistBins=40;

for vaIdx = 1:length(varargin)
    currArg = varargin{vaIdx};
    notLastArg = vaIdx~=length(varargin);
    switch(currArg)
        case {'mean','Mean'}
            meanFlag=1;
        case {'median','Median'}
            medianFlag=1;
        case {'max','Max'}
            maxFlag=1;
        case {'min','Min'}
            minFlag=1;
        case {'hist','Hist'}
            histFlag=1;
            if notLastArg && ~ischar(varargin{vaIdx+1})
                histBins=varargin{vaIdx+1};
            end
        case {'varHist','VarHist','varhist','Varhist'}
            varHistFlag=1;
            if notLastArg && ~ischar(varargin{vaIdx+1})
                varHistBins=varargin{vaIdx+1};
            end
        case {'CI','Ci','ci'}
            if notLastArg && ~ischar(varargin{vaIdx+1})
                confInterval=varargin{vaIdx+1};
            else
                disp('simulateMC(): Value is missing from ''CI'' name-value pair.')
            end
        otherwise
            if ischar(currArg)
                disp("simulateMC(): The string '"+currArg+"' in the function parameters was not recognized as a valid option.")
            end
    end
end

%% Generate & Store Samples From the Specified Distributions
distIdx=0;
sampIdx=0;
while(distIdx < size(dists,1))
    distIdx=distIdx+1;
    sampIdx=sampIdx+1;
    distType = dists{distIdx}{1};
    switch(distType)
        case {'Bootstrap','bootstrap'}
            %Bootstrap a sample of size N from user provided data.
            %Syntax is:  {'Bootstrap',datavec}
            data=dists{distIdx}{2};
            [rown,coln]=size(data);
            if coln~=1
                data=transpose(data);
            end
            sampledVars{sampIdx} =  randsample(data,N,true);
        case {'BootstrapMean','bootstrapMean','Bootstrapmean','bootstrapmean'}
            %Boostrap an approximate distribution of the sample mean of
            %user provided data.
            %Syntax is:  {'BootstrapMean',datavec}
            data=dists{distIdx}{2};
            [rown,coln]=size(data);
            if coln~=1
                data=transpose(data);
            end
            means=zeros(N,1);
            for i=1:N
                means(i)=mean(randsample(data,length(data),true));
            end
            sampledVars{sampIdx} = means;
        
        case {'Corr','corr'}
            %Use an iteratively optimized gaussian copula approach to 
            %generate correlated samples from the specified distributions.
            %Syntax is:  {'Corr',distParamLists,PearsonCorrelationMatrix,tolerance}
            
            %The tolerance parameter is optional. The default tolerance value 
            %specifies a maximum error of 1% in the correlation matrix of 
            %the generated samples.
            corrDistParams=dists{distIdx}{2};
            numDists = length(corrDistParams);
            for i=1:numDists
                corrDists(i)=getDistribution(corrDistParams{i});
            end
            
            corrMatGoal = dists{distIdx}{3};
            if length(dists{distIdx})>3
                tol = dists{distIdx}{4};
            else
                tol = 0.01;
            end
            
            %Iteratively optimize the copula's correlation matrix until the 
            %correlation matrix of the resulting samples is within tolerance.
            maxErrorAbs=1;
            corrMat=corrMatGoal;
            means = zeros(1,numDists);
            alpha = 0.1;    %iterative update ratio, adjusting alpha may improve convergence
            seed = clock; seed=seed(5)*seed(6)^2*rem(seed(5),seed(6));
            maxIters=10000; count=0;
            while maxErrorAbs>tol
                rng(seed)  %using the same rng seed each iteration can aid convergence
                copula = mvnrnd(means,corrMat,N);
                copula = normcdf(copula);
                samps = zeros(N,numDists);
                for i=1:numDists
                    samps(:,i) = icdf(corrDists(i),copula(:,i));
                end
                corrErrorMat = (corrMatGoal-corr(samps))./corrMatGoal;
                maxErrorAbs=max(max(abs(corrErrorMat)));
                [rownum,colnum]=find(abs(corrErrorMat)==maxErrorAbs);
                maxError = corrErrorMat(rownum(1),colnum(1));
                corrMat(rownum(1),colnum(1))=(1+(maxError*alpha))*(corrMat(rownum(1),colnum(1)));
                corrMat(colnum(1),rownum(1))=(1+(maxError*alpha))*(corrMat(colnum(1),rownum(1)));
                %TODO: Implement a check to ensure corrMat can't be
                %repeatedly updated in the wrong direction. 
                count=count+1;
                if count>maxIters
                   disp("simulateMC(): The correlated samples were not obtained within tolerance in "+maxIters+" iterations.")
                   return
                end
            end
            for i=1:numDists
                sampledVars{sampIdx+i-1} = samps(:,i);
            end
            sampIdx = sampIdx + numDists - 1;  
            
        case {'Custom','custom'}
            %Let the user provide the data sample
            %Syntax is: {'Custom',datavec}
            data = dists{distIdx}{2};
            [rown,coln]=size(data);
            if coln~=1
                data=transpose(data);
            end
            [rown,coln]=size(data);
            if rown~=N
                disp("simulateMC(): The length of the custom sample must equal the PropUncertMC parameter 'N'.")
                return
            end
            sampledVars{sampIdx} = data;
        case {'Fit','fit'}
            %Fits the specified type of distribution to user provided data
            %and generates a sample of length N from the fitted distribution.
            %Syntax is: {'Fit','DistributionType',data}
            fittedDist= fitdist(dists{distIdx}{3},dists{distIdx}{2});
            sampledVars{sampIdx} = getSamples(fittedDist,N);
        case {'Uniform','uniform'}
            %Enables uniform samples to be generated without the need for 
            %the Statistics and Machine Learning Toolbox.
            %Syntax is: {'Uniform',upperbound,lowerbound}
            sampledVars{sampIdx} = dists{distIdx}{2} + (dists{distIdx}{3}-dists{distIdx}{2})*rand(N,1);
        otherwise
            %Sample a distribution specified by the distribution params
            %   Syntax is: {'DistributionName',param_1, ... ,param_n}
            %One can also specify a truncated distribution like so:
            %   {'DistributionType',param_1, ... ,param_n,'trunc',lowerCutOff,upperCutOff}
            
            if ~ischar(distType)
                %If no distribution type is specified then assume normal
                dists{distIdx}=horzcat('Normal',dists{distIdx});
                distType = 'Normal';
            end
            if strcmp(distType,'Normal') || strcmp(distType,'normal')
                %Enables normal samples to be generated without the need
                %for the Statistics and ML toolbox.
                truncIdx=getTruncIdx(dists{distIdx});
                if truncIdx<1
                    sampledVars{sampIdx} = dists{distIdx}{2} + dists{distIdx}{3}*randn(N,1);
                else
                    sampledVars{sampIdx} = getTruncatedNormals(dists{distIdx}{2},dists{distIdx}{3},dists{distIdx}{truncIdx+1},dists{distIdx}{truncIdx+2},N);
                end
            else
                
            %Handles the sampling for distributions besides normal and
            %uniform.
            sampledVars{sampIdx} = getSamples(getDistribution(dists{distIdx}),N);
            end
    end
end


%% Evaluate Function & Get Descriptors
MCsamples=cell2mat(sampledVars);
MCfuncVals=func(sampledVars{:});
stdev = std(MCfuncVals);

%Find smallest interval containing 'confInterval' of the values
sorted = sort(MCfuncVals);
range = ceil(confInterval*length(sorted));
smallestInterval = inf;
bottom=-inf; top=inf;
for idxLow = 1:length(sorted)
   idxHigh = idxLow+range;
   if idxHigh > length(sorted)
       break
   end
   low = sorted(idxLow);
   high = sorted(idxHigh);
   if (high-low)<smallestInterval
       bottom = low; top = high;
       smallestInterval = high-low;
   end
end
CI = [bottom top];
if (1-confInterval)*length(sorted) < 1000
   recomm = floor(1000/(1-confInterval));
   disp("simulateMC(): Recommended to use at least "+recomm+" samples for a CI interval this large")
end

% Set funcVal as specified by the optional arguments
if meanFlag
    funcVal = mean(MCfuncVals);
elseif medianFlag
    funcVal = median(MCfuncVals);
elseif maxFlag
    funcVal = max(MCfuncVals);
elseif minFlag
    funcVal=min(MCfuncVals);
else
    funcVal = mean(CI);
end

%% Plot histograms as specified by the optional arguments
if varHistFlag
    colN=ceil(sqrt(sampIdx));
    rowN = ceil(sampIdx/colN);
    figure
    for pIdx = 1:sampIdx
        subplot(rowN,colN,pIdx);
        hist(cell2mat(sampledVars(:,pIdx)),varHistBins);
        title("Var "+pIdx)
    end
    sgtitle('Histogram Of MC Samples For Each Variable')
end

if histFlag
    figure
    hist(MCfuncVals,histBins)
    title('Histogram Of MC Function Values')
end


%% Helper Functions
    function distrib = getDistribution(input)        
        if ~ischar(input{1})
            %assume normal if not specified
            input=horzcat('Normal',input);
        end
        
        %truncate distribution if neccesary
        truncInd=getTruncIdx(input);        
        if truncInd<1
            distrib = makedist(input{:});
        else
            distrib = makedist(input{1:truncInd-1});
            distrib = truncate(distrib, input{truncInd+1},input{truncInd+2});
        end
    end
        
    
    function tIdx = getTruncIdx(list)
        tIdx=0;
        truncWords = {'t' 'T' 'Trunc' 'Truncate' 'trunc' 'truncate'};
        for idx = 1:length(list)
            if any(strcmp(list{idx},truncWords))
                tIdx = idx;
            end
        end
    end


    function samples = getSamples(distrib,num)
        samples = random(distrib,num,1);
    end


    function tNorm = getTruncatedNormals(Mu,Sigma,lowerBound,upperBound,n)
            tNorm =  Mu + Sigma*randn(n,1);
            tooHigh=numel(tNorm(tNorm>upperBound)); 
            tooLow=numel(tNorm(tNorm<lowerBound));
            while(tooHigh>0 || tooLow>0)
                if(tooHigh>0)
                    tNorm(tNorm>upperBound) = Mu + Sigma*randn(numel(tNorm(tNorm>upperBound)),1);
                end
                if(tooLow>0)
                    tNorm(tNorm<lowerBound) = Mu + Sigma*randn(numel(tNorm(tNorm<lowerBound)),1);
                end
                tooHigh=numel(tNorm(tNorm>upperBound)); 
                tooLow=numel(tNorm(tNorm<lowerBound));
            end
    end
end
    