function barrierRankYalmip()

if ~checkDependency('yalmip')
  error('Yalmip not found');
end

% Variables/Indeterminates
sdpvar x1 x2
x = [x1;x2];

% Dynamics
f = [x2; -x1 + x1^3 - x2];

% Description of initial condition set
g_X0 = 0.5^2 - (x1 - 1.5)^2 - x2^2; % Should be in terms of x1 and x2

% Description of unsafe set
g_Xu = 0.4^2 - (x1+1)^2 - (x2+1)^2; % Should be in terms of x1 and x2

% Multipliers
dL1=2;
[L1,cL1,mL1]=polynomial(x,dL1);

dL2=2;
[L2,cL2,mL2]=polynomial(x,dL2);

dL=2;
[L,cL,mL]=polynomial(x,dL);

% Barrier function
[B,cB,mB] = polynomial(x,4);

% Compute time derivative of B
Bdot = jacobian(B,x)*f;

constr = -Bdot + L*B; % Bdot < 0 when B = 0
[cconstr, mconstr] = coefficients(constr,x);

w = [cB;cL];
W = sdpvar(length(w));
ww = w*w';
%sdisplay(cVdot)
for i = 1:length(W)
  for j = 1:length(W)
    aux = getbasematrix(cconstr,getvariables(ww(i,j)));
    cconstr = cconstr-aux*ww(i,j)+aux*W(i,j); %remember you are rewriting.
  end
end
constr=cconstr.'*mconstr; %rewriting constr.


M=[W w;w' 1];

FF=[sos(constr), M>0, ... % bilinear constraint
  sos(-B - L1*g_X0 - 1e-2), ... % initial condition set constraint
  sos(B - L2*g_Xu), ... % unsafe set constraint
  sos(L1), sos(L2), ... constrain multipliers to be sos
  ];


sol_w=zeros(length(w),1);
for k=1:20 %this loop improves on the standard trace heuristic for rank.
  solvesos(FF,trace(W)-2*sol_w'*w,[],[cB;cL;cL1;cL2;vec(W)])
  %solvesos(FF,trace(W)-2*sol_w'*w,[],[cB;cL;vec(W)])
  sol_w=double(w);
%   eig(double(W))%just so we keep track of it.
end
sol_w=double(w)
sol_W=double(W);
eig_sol_W=eig(sol_W) %if rank 1, we are good, W must be w*w'.
residual = sol_w*sol_w'-sol_W


sol_B=double(cB).'*mB;
%sdisplay(clean(sol_B,1e-7))
sdisplay(sol_B)


% Plot stuff
figure
[X,Xd] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));
% initial condition set
gs_X0 = subs(g_X0,x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_X0),100,100),[0 0],'g','LineWidth',3); % initial condition set
hold on
% unsafe set
gs_Xu = subs(g_Xf,x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_Xu),100,100),[0 0],'r','LineWidth',3) % unsafe set
% 0-level set of B
gs_B = subs(sol.eval(B),x,[X(:),Xd(:)]');
contour(X,Xd,reshape(double(gs_B),100,100),[0 0],'b','LineWidth',3) % barrier zero level set

% (scaled) vector field
[X,Xd] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
x1dots = reshape(double(msubs(f(1),x,[X(:),Xd(:)]')),50,50);
x2dots = reshape(double(msubs(f(2),x,[X(:),Xd(:)]')),50,50);
x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
quiver(X,Xd,x1dots,x2dots,'AutoScale','off','Color','k');

% Title
title('Barrier functions')
xlabel('x_1');
ylabel('x_2');
legend('X_0','X_u','0 level-set of B(x)');

end
