function [X,IM,D]=GLMDAUG(X0,beta,m,nruns,model,family,link, RP,NRP)
%
% DGLMaug Augments a Local D-optimal Experimental Design for 
% multivariate Generalized Linear Models (GLM) to a total of nruns
% (original + new design points)
% using sequential Monte Carlo 
% and Federov's Exchange Algorithm as implemented in MATLAB statistical toolbox (required)
%
% Syntax:
% X=GLMDAUG(X0,beta,m)
% [X,IM,D]=GLMDAUG(X0,beta,m,nruns,model,family,link, RP,NRP)
%
% Description:
% X=GLMDAUG(X0,beta) is the minimal form, which adds one point to a local D-optimal design
%                                   for a 'binomial' First-Order model with the logit link
%                                   with number of variables, m, determined
%                                   by X0 - The fixed design points
%                                   All variabels assumed to be constrained to the region [-1,1]
%                                   beta - a vector containing the  coefficient values
%                                   X is a matrix containing the output m  variables
%                                        values at the suggested experimental design
%
% Optional input variables are:
% m: number of variables (important if model is not linear!)
% nruns: number of runs in the final design; a scalar. If [], default=2^m
% model: a model in the form legitimate for the function x2fx
%                 options are 'linear' (default), 'interaction',
%                 'quadratic', 'purequadratic', or a MATRIX of model terms
%                  A matrix example: model=[0 0 0 ; 1 0 0; 0 1 0 ; 0 0 1] is 
%                      equivalent to model='linear', for m=3
% family: implemented families are: 'binomial' (default), and 'poisson'
% link: 'logit' (default for binomial), 'probit', 'cll', 'log' (default for Poisson)
% RP: Required Precision (default 10^-2) of X matrix elements
% NRP: Number of random points around each base (best guess) of X from the
%           previous iteration step (default = 50)
%
% Values not supplied must be replaced by []
% for example X=DGLMaug(X0,[0 2 2 0],2,[],'interaction',[],'CLL',10^-4)
%
% [X,IM]=GLMDAUG(...)   returns, in the matrix IM, the information matrix of [X0;X]
% [X, IM,D]=GLMDAUG(...)   returns in the scalar D the D-criteria value for
%                                       the model; that is the p-th root of the determinant of 
%                                       the information matrix IM, divided  by the number of runs (nruns).
%                                       p is the length of the vector beta
%                                       (the number of model coefficients)
%
% This script was written by Hovav Dror, Tel-Aviv University, 2006
%
% If you make changes to the script, or if you implement it
% to a (to be) published problem, please inform us:
% dms@post.tau.ac.il ; hovavdror@gmail.com
%
% further details, including papers describing the algorithm, can be
% found at: http://www.math.tau.ac.il/~dms/GLM_Design
 

t=clock;
 
if nargin<9, NRP=[];
    if nargin<8, RP=[];
        if nargin<7, link=[];
            if nargin<6, family=[];
                if nargin<5, model=[];
                    if nargin<4, nruns=size(X0,1)+1;
                        if nargin<3, m=size(X0,2)-1;
                    end, end, end, end, end, end, end
 
if nargin==0
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf('[X,IM,D]=DGLMaug(X0,beta,m,nruns,model,family,link, RP,NRP) \n');
        fprintf('Type "help DGLMaug" or see the function''s code\n');
        fprintf('for help and details \n');
%        fprintf('use: help DGLM \n');
        fprintf('----------------------------------------------------------------------------------\n');
elseif nargin<2
    fprintf('Error: Not enough arguments.  \n');
else % ----------- if there are enough parameters, we can begin
    
    % ---------- substitute defaults ------------------
    if isempty(nruns), nruns=2^m;end;
    NeededPoints=nruns-size(X0,1); % Number of Needed points to augment
    if isempty(model), model='linear'; end;  
    if isempty(family), family='binomial'; end;
    if isempty(link)
        if strcmp(family,'poisson'), link='log';
        else link='logit'; end, end
    if isempty(RP), RP=10^-2; end;
    if isempty(NRP), NRP=50; end;
        
    if isnumeric(model) % Decide what is p
            p=length(model);
     else
            p=length(x2fx(ones(1,m),model));
     end; % if isnumeric(model)
     if size(beta,1)==1, beta=beta'; end;
     if p > length(beta)
            fprintf('---------------------------------------------------------------------------------------------------\n');
            fprintf('Error: Size of Beta is shorter than the requested model ! \n');
            fprintf('---------------------------------------------------------------------------------------------------\n');
     else
            if p<length(beta)
                beta=beta(1:p);
                fprintf('Attention: beta is too long. Truncating it. \n');
           end;
    end;  % if p~=length(beta)
    if nruns<p
        fprintf('Number of runs too small, increasing it');
        nruns=p;
    end;
 
    
    % ----------- More initialization ---------------
        if size(X0,1)>0
            F0=x2fx(X0,model); % Turn X0 into a design matrix (without the weights, yet)
            F0=diag( sw(F0,beta,family,link) ) * F0; % Putting the (square root of ) weights
        else
            F0=[];
        end;

    
    SearchR=1; % Initial Search Radius
    
    X=rand(NeededPoints,m)*2-1; % Initial random design
    X=sortrows(X);
    if (NeededPoints*NRP) < 1000  % ensure initial grid is big enough
        CS=rand(1000-NeededPoints*NRP,m)*2-1;    % CS: Candidate Set
    else
        CS=[];
    end;
    % ---------------- Begin Search Loop------------------------
    while SearchR>RP % while the Required Precision not achieved
        CS=[CS; X];
        for i=1:NeededPoints
            CS=[CS; ones(NRP,1)*X(i,:) + randn(NRP,m)*SearchR];% Add random points around X
        end;
        CS(CS < -1) = -1;% ensure the variables limits are [-1,1]
        CS(CS > 1) = 1 ;
        
        CSF=x2fx(CS,model); % Turn into a design matrix (without the weights, yet)
        WCSF=diag( sw(CSF,beta,family,link) ) * CSF; % Putting the (square root of ) weights
                
        % -------- Using the "regular" D-optimal algorithm ---------
        rows=canddaugm(F0,WCSF,nruns,WCSF(1:NeededPoints,:)); %Choose the D-optimal nruns from C. return their row numbers.
        % WCSF is the Candidate Set, its first nruns rows are the initial design
 
        % ----------- Updating the Search Radius
        XN=sortrows(CS(rows,:));
        G=XN - X;
        DiffG=max(max(abs(G)));
        SearchR1=min( [DiffG   SearchR]);
        if (SearchR1==SearchR)  % if no improvement in Search radius, check if new design really better
            XN2=x2fx([X0;XN],model);
            X2=x2fx([X0;X],model);
            XN2=diag( sw(XN2,beta,family,link) )*XN2;
            X2=diag( sw(X2,beta,family,link) )*X2;
            if (det(XN2'*XN2) / det(X2'*X2)) < 1 % if no improvement
                XN=X; % retain old design
            end;
        end;
        SearchR=max([SearchR1 0.3*SearchR]);% Ensuring no RMradius=0 or too small change in RMradius
        CS=[]; % Empty the Candidate Set
       
        X=XN; % Updating the best design found so far
        
        if etime(clock,t)>10 % every 10 seconds show a sign of life with accuracy report
            t=clock;
            fprintf('Current Search Radius %g \n', SearchR);
        end;
        
    end; % While SearchR>RP
    
end; % if nargin==0
 
% -------------------- Calculate Results for IM and D -----------------------------------
if nargout>1
    FX=x2fx([X0 ; X],model);
    IM=diag( sw(FX,beta,family,link) ) * FX;
    IM = IM' * IM;
    if nargout>2
        D=det(IM)^(1/p) / nruns; 
    end;
end; % if nargout>1
 
    
    % --------------- nested functions ----------------------------
    function y=sw(F,b,family,link) % Decision on (Sqrt of) Weights
        if strcmp(family, 'binomial')
            if strcmp(link,'logit')
                y=bsw(F,b);
            elseif strcmp(link,'CLL')  | strcmp(link,'comploglog');
                y=cllsw(F,b);
            elseif strcmp(link,'probit')
                y=probitsw(F,b);
            else
                fprintf('----------------------------------------------------\n');
                fprintf('Error: Unkown link function \n');
                fprintf('----------------------------------------------------\n');
           end; % if link
        elseif strcmp(family,'poisson')
            y=sqrt(exf(F,b)); % assuming log link automaticaly
        else
                fprintf('----------------------------------------------------\n');
                fprintf('Error: Unkown family function \n');
                fprintf('----------------------------------------------------\n');          
        end;
    end % sw
    
    % ---------------  Additional Support Functions -----------------------
    function y=bsw(F,b)
        ex=exf(F,b);
        y= sqrt (ex ./ (1+ex).^2 );
    end
    function y=cllsw(F,b)
            ex=exf(F,b);
            y=sqrt(ex.*ex.*exp(-ex) ./ (1.-exp(-ex)) );
    end
    function y=probitsw(F,b)
            p=normcdf((F*b)')';
            y= sqrt( 1./(p.*(1-p)) .* exp(-0.5*(F*b).^2) ./ (2*pi).^0.5 );
    end
        
    function y=exf(F,b)
        y=exp(F*b);
    end
 % ------------------------------------------ end of support nested functions ------------------------------
 
 
end % DGLM function
 

 