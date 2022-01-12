function IM=InfoMtrxGLM(F,beta,family,link)
% A function calculating the Infomation Matrix
% for a regression matrix F ( F=x2fx(X,model))
% given a GLM model with specified family and link
if size(beta,2)>1, beta=beta'; end;
if strcmp(family,'binomial')
    if strcmp(link,'logit')
        EX=exp(F*beta);
        W=EX ./ (1+EX).^2 ; 
     elseif strcmp(link,'cll') | strcmp(link,'comploglog')
        EX=exp(F*beta);
        W=EX.*EX.*exp(-EX) ./ (1.-exp(-EX)) ;
    elseif strcmp(link,'probit')
        p=normcdf((F*beta)')';
        W=  1./(p.*(1-p)) .* (exp(-0.5*(F*beta).^2) ./ (2*pi).^0.5).^2 ;
    else
        fprintf('Error! Unknown link');
    end;
elseif strcmp(family,'poisson')
    if strcmp(link,'log')
        W=exp(F*beta);
    else
        fprintf('Error! Unknown link');
    end;
else
    fprintf('Error: Unknown family');
end;
        
IM=F'*diag(W)*F;