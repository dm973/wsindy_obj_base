%%% how to construct a variety of terms
%% function composition
t1 = term('ftag',[1 0 3]);
disp(t1.get_str)

t2 = term('fHandle',@(x)atan(x));
disp(t2.get_str)

t3 = compterm(t2,t1);
disp(t3.get_str)

%% vector function composition: 1 layer
t1 = library('tags',{[1 2 0 3],[0 1 1 0],[0 0 1 -1]});
cellfun(@(t)disp(t.get_str),t1.terms)

t2 = term('fHandle',@(x,y,z)exp(-x.^2-y.^2).*z);
disp(t2.get_str)

t3 = compvec(t2,t1);

N = 1000;
X = sin(linspace(0,2*pi,N)'*randi(10,1,4));
Y = t3.evalterm(X);

plot(Y)

%% vector function composition: 2 layers
t1 = library('tags',{[1 2 0 3],[0 1 1 0],[0 0 1 -1]});
cellfun(@(t)disp(t.get_str),t1.terms)

t2 = library('tags',{[1 1 0],[0 1 1]});
cellfun(@(t)disp(t.get_str),t2.terms)

t3 = library('terms',cellfun(@(t)compvec(t,t1),t2.terms,'un',0));

t4 = term('fHandle',@(x,y)exp(-x.^2-y.^2));

t5 = compvec(t4,t3);

N = 1000;
X = sin(linspace(0,2*pi,N)'*randi(10,1,4));
Y = t5.evalterm(X);

plot(Y)

%% vector function composition: 2 layers - single input

t1 = library('tags',-1i*[1:6]');
t2 = library('tags',get_tags(0:2,[],6));
t3 = library('terms',cellfun(@(t)compvec(t,t1),t2.terms,'un',0));
t4 = term('fHandle',@(varargin)exp(-sum([varargin{:}].^2,2)),'nstates',length(t3.terms));
t5 = compvec(t4,t3);

N = 1000;
x = linspace(0,2*pi,N)';
X = sin(x*randi(20));
Y = t5.evalterm(X);

plot(x,X,x,Y)

%% gradient computations

% t1 = term('ftag',[1 2],'linOp',[2 0]);
% t2 = term('ftag',3);
% t3 = compterm(t2,t1);

t1 = term('fHandle',@(u,v)sin(u-v),'linOp',[1 0]);
t2 = term('fHandle',@(u,v)v.*2,'linOp',[1 0]);
t3 = prodterm(t1,t2);

t3.get_grads;
t3.get_str
t3.gradterms(1).get_str
t3.gradterms(2).get_str

Nx = 80; Nt = 100;
x = linspace(-1,2,Nx)';
t = linspace(0,10,Nt);
u = sin(x).*t;
v = atan(x-t);

Uobj = wsindy_data({u,v},{x,t});

Y = t3.evalterm(Uobj);
Y1 = t3.gradterms(1).evalterm(Uobj);
Y2 = t3.gradterms(2).evalterm(Uobj);

subplot(321)
imagesc(u)
title('u')
colorbar
subplot(322)
imagesc(v)
title('v')
colorbar
subplot(323)
imagesc(Y1)
title('D_u Y')
colorbar
subplot(324)
imagesc(Y2)
title('D_v Y')
colorbar
subplot(325)
imagesc(Y)
title('Y')
colorbar

%%


%{prodterm(term('ftag',1),diffOp([1 0],'meth','wffd'))};
% custom_terms = {compterm(term('ftag',2),diffOp([1 0],'meth','wffd'))};
