% Detailed examples for the function simulateMC()

% A list of supported distributions and their parameters can be found here:
%      https://www.mathworks.com/help/stats/supported-distributions.html
%Joe Klebba, 3/2021

%EXAMPLE 1
%Combine samples from X and Y through the function f(X,Y)=5*x+y. 
%X is a normal distribution with params [mean=5,sigma=1]. 
%Y is a uniform distribution with params [lowerbound=6,upperbound=7].
clear, clc
disp('Example 1')
f=@(X,Y)5.*X+Y;
dists = {{'Normal',5,1};
         {'Uniform',6,7}};
n=100000;
[funcVal,stdev,CI] = simulateMC(f,dists,n)

%Call the function again but this time provide the 'varHist' name-value 
%argument in order to plot a 65 bin histogram of each variable's sample
%Also provide the 'hist' argument in order to plot the resulting 
%distribution of f(X,Y,Z). Use the default number of bins.
[funcVal,stdev,CI] = simulateMC(f,dists,n,'varHist',65,'hist');
return


%EXAMPLE 2
%Combine samples from X, Y, and Z through the function f(X,Y,Z)=X^2/(Y*Z).
%X is a truncated normal distribution with params [mean=4,sigma=0.2,lowerb=3.5,upperB=4.5]. 
%Y is a triangular distribution that is correlated with Z. Params are [lowerB=3.3,peak=3.5,upperB=3.7].
%Z is a lognormal distribution that is correlated with Y. Params are [mean=1,sigma=0.05].
%The correlation betwen Y and Z is specified by the Pearson correlation matrix 'corrMat'.

disp('Example 2')
f=@(X,Y,Z)X.^2./(Y.*Z);
corrMat = [1 0.5; 0.5 1];
dists = {{'Normal',4,0.2,'trunc',3.5,4.5};
         {'Corr',{{'Triangular',3.3,3.5,3.7};
                  {'Lognormal',1,0.05}},corrMat}};
n=100000;
%Use 'CI' argument to specify a pseudo confidence interval of 0.95.
%Use 'mean' argument to return the mean of f(X,Y,Z) as funcVal instead 
% of using the default.
[val,stdev,CI,MCfuncVals,MCsamples] = simulateMC(f,dists,n,'CI',0.95,'mean');

%Also find the median value of the resulting distribution.
med = median(MCfuncVals)

%Show that the samples of Y and Z have approximately the right correlation (0.5).
checkCorr = corrcoef(MCsamples(:,2),MCsamples(:,3))
clear MCfuncVals MCsamples %clear these to save memory 
return


%EXAMPLE 3
%Combine samples from W,X,Y,Z through the function f.
%For W a sample will be bootstrapped from provided data.
%For X a normal distribution will be fit to provided data and samples will
%be drawn from the fitted distribution.
%For Y a custom data sample will be passed (it must be of length N).
%Z will be a Weibull distribution with params [scale=3,shape=2].
disp('Example 3')
n=100000;
Wdata = [5 5.1 5.1 5.2 5.3 5.3 5.2 5.1 5.2 5.2 5.3 5.4 4.9 5.3];
Xdata = normrnd(20,1,50,1);
Zdata = unifrnd(10,10.5,n,1);

f=@(W,X,Y,Z)W.*X+Y+5.*Z;
dists = {{'Bootstrap',Wdata}
         {'Fit','Normal',Xdata}
         {'Custom',Zdata}
         {'Weibull',3,2}};
funcVal = simulateMC(f,dists,n);
return
%Note: The methods 'Bootstrap','Fit','Custom', and 'BootstrapMean' cannot be
%utilized by simulateMC() for generating correlated samples.


%EXAMPLE 4
%Combine samples from six different distributions where a group of 3
%distributions are correlated and another group of 2 are correlated. 
%Some of the distributions are truncated.

%Specify a maximum correlation error of 0.1% for the correlated
%distributions. (The default max error is 1%).

%Plot histograms of all the samples and the function outputs.
disp('Example 4')
tol=0.001;
corrMat1 = [1 0.3 -0.7; 
            0.3 1 -0.5; 
            -0.7 -0.5 1];
corrMat2 = [1 0.2; 0.2 1];
dists =  {{'Corr',{{'normal',5,0.1,'trunc',4.4,5.2};
                   {'uniform',4,4.001};
                   {'uniform',3,3.002}},corrMat1,tol};
          {'triangular',10.1,10.25,10.27};
          {'Corr', {{'normal',4,0.2};
                    {'triangular',7,8,9,'trunc',7.5,inf}},corrMat2,tol}};
f=@(x,y,z,q,p,r)x+y.*z+q+p-r;
funcVal=simulateMC(f,dists,n,'varHist',70,'hist',70);

