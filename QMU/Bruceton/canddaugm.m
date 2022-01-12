function rowlist = canddaugm(F0,fxcand,nruns,xinit)

% This function augments a D-optimal design using the row exchange algorithm
% It is based on the function CANDEXCH from the statistical toolbox
% modifying it so it can be used to augment designs

ReadySize=size(F0,1);
AugmSize=nruns-ReadySize;
rowlist = zeros(AugmSize,1);
maxiter = 10;
iter = 0;
madeswitch = 1;
dcutoff = 1 + sqrt(eps(class(fxcand)));

okargs = {'maxiter' 'init' 'display'};
defaults = {10 [] 'on'};
maxiter=50;
dodisp='off';

quiet = isequal(dodisp,'off');
if isempty(xinit)
   X = fxcand(unidrnd(size(fxcand,1),AugmSize,1),:);
   X=[ X ; F0];
else
   X = xinit;
   X=[ X ; F0];
end
F = [];

[Q,R]=qr(X,0);

% Adjust starting design if it is rank deficient, because the algorithm
% will not proceed otherwise.
p = size(R,2);
if rank(R)<p
   warning('stats:candexch:BadStartingDesign',...
           'Starting design is rank deficient');
   diagr = abs(diag(R));
   smallval = sqrt(eps(class(diagr))) * max(diagr);
   t = (diagr < smallval);
   if any(t)
      tind = (1:p+1:p^2);
      R(tind(t)) = smallval;
   end
end
logdetX = 2*sum(log(abs(diag(R))));

% Create Iteration Counter Figure.
if ~quiet
   screen = get(0,'ScreenSize');
   f = figure('Units','Pixels','Position',[25 screen(4)-150 300 60],...
              'MenuBar','none','Name','D-Optimal Design Generation',...
              'NumberTitle','off');
   ax = gca;
   set(ax,'Visible','off');
   t=text(0,0.7,'Row exchange iteration counter.');
   set(t,'FontName','Geneva','FontSize',12);
   itertxt = 'Iteration %4i';
   s=sprintf(itertxt,1);
   set(f,'CurrentAxes',ax);
   h=text(0,0.2,s);
   set(h,'FontName','Geneva','FontSize',12);
end

while madeswitch > 0 & iter < maxiter
   madeswitch = 0;
   iter = iter + 1;

   % Update iteration counter.
   if ~quiet & ishandle(h)
      s=sprintf(itertxt,iter);
      set(h,'String',s);
      drawnow;
   end
  
   for row = 1:AugmSize
      % Compute determinant factor over whole candidate set if not done yet
      if isempty(F)
         F = fxcand/R;
         dxcand = sum(F.*F, 2);
      end
      
      E = X(row,:)/R;
      dxold = E*E';
      dxswap  = F*E';
      
      dd = (1+dxcand) * (1-dxold) + dxswap.^2;
     
      % Find the maximum change in the determinant.
      [d,idx] = max(dd);
     
      % Switch rows if the maximum change is greater than 1.
      if (d > dcutoff) | (rowlist(row) == 0)
         madeswitch = 1;
         logdetX = log(d) + logdetX;
         rowlist(row) = idx;
         X(row,:) = fxcand(idx,:);
         [Q,R] = qr(X,0);
         F = [];  % needs re-computing using new R
      end
   end
end
if ~quiet & ishandle(f), close(f); end
   
