clear all;
%%
sr = load('sim/base/sub01_x3932160x_x1835008x_base.mat');sr = sr.obj;
p = sr.fit(end).p;

% p.v = 2;
p.st = 0.1;
p.sv = 2.5;
[pdf_cw_bru,pdf_cr_bru,rt_bru,cdf_dow,cdf_ups] = sr.ddm_bru(p,sr.fit(end).s.dt,sr.fit(end).s.T,25e3);
%%

%%
p2.a = p.a*2;
p2.v = p.v;
p2.t = p.t;
p2.st = p.st*2;

p2.sz = 0;
p2.sv = p.sv;
p2.z = 0.5;
%%
err = 1e-5;
rfpdf = @(x) full_pdf(x,p2.v,p2.sv,p2.a,p2.z ,p2.sz ,p2.t,p2.st,err);
op = prob_ub(p2.v,p2.a,p2.z);
rt_r = +[0:1e-3:2];
rt_w = -[0:1e-3:2];

% clf;hold on;
%%
[pdf_cw_trm,pdf_cr_trm,rt_trm] = sr.ddm_pdf(p,sr.fit(end).s.dt,sr.fit(end).s.T,250);

%%
clf;hold on;
% subplot(2,1,1);cla;hold on;
% subplot(2,1,2);cla;hold on;
% full_pdf(0.3,p.v,0,p.a,0 ,0 ,p.t,p.st,err)
p_pdf = arrayfun(@(rt) rfpdf(rt),rt_r);
% subplot(2,1,1);
% plot(rt_r,p_pdf*(1-op));
% subplot(2,1,2);
plot(rt_r,p_pdf,'b');

p_pdf = arrayfun(@(rt) rfpdf(rt),rt_w);
% subplot(2,1,1);
% plot(rt_r,p_pdf*(op),'b--');
% subplot(2,1,2);
plot(rt_r,p_pdf,'b--');

plot(rt_bru,pdf_cr_bru,'r');
plot(rt_bru,pdf_cw_bru,'r');

plot(rt_bru,pdf_cw_trm,'g');
plot(rt_bru,pdf_cr_trm,'g');
xlim([0 2])
%%
function p = ftt_01w(tt,w,err)
if (pi*tt*err)<1 % if error threshold is set low enough
    kl=sqrt(-2*log(pi*tt*err)/(pi^2*tt)); % bound
    kl=max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met
else % if error threshold set too high
    kl=1/(pi*sqrt(tt)); % set to boundary condition
end
% calculate number of terms needed for small t
if (2*sqrt(2*pi*tt)*err)<1 % if error threshold is set low enough
    ks=2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound
    ks=max(ks,sqrt(tt)+1); % ensure boundary conditions are met
else % if error threshold was set too high
    ks=2; % minimal kappa for that case
end

% compute f(tt|0,1,w)
p=0; %initialize density
if ks<kl % if small t is better (i.e., lambda<0)
    K=ceil(ks); % round to smallest integer meeting error
    vec_k = -floor((K-1)/2):ceil((K-1)/2);
    for k = vec_k
        p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
    end
    p=p/sqrt(2*pi*(tt^3)); % add con_stant term
    
else % if large t is better...
    K= ceil(kl); % round to smallest integer meeting error
    
    for k=1:K
        p=p+k*exp(-((k^2))*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
    end
    p=p*pi; % add con_stant term
end
end

function op = prob_ub(v,a,z)
%"""Probability of hitting upper boundary."""
op = (exp(-2*a*z*v) - 1) / (exp(-2*a*v) - 1);
end


function op = pdf(x,v,a,w,err)
%"""Compute the likelihood of the drift diffusion model f(t|v,a,z) using the method
%and implementation of Navarro & Fuss, 2009.
%"""
if x <= 0
    op = 0;
else
    
    tt = x/(a^2);% # use normalized time
    p = ftt_01w(tt, w, err);% #get f(t|0,1,w)
    
    %# convert to f(t|v,a,w)
    op =  p*exp(-v*a*w -((v^2))*x/2.)/((a^2));
end
end

function op = pdf_sv(x,v,sv,a,z,err)
%"""Compute the likelihood of the drift diffusion model f(t|v,a,z,sv) using the method
%and implementation of Navarro & Fuss, 2009.
%sv is the std of the drift rate
%"""
if x <= 0
    op = 0;
elseif sv==0
    op =  pdf(x, v, a, z, err);
else
    tt = x/((a^2)); %# use normalized time
    p  = ftt_01w(tt, z, err); %#get f(t|0,1,w)
    
    %# convert to f(t|v,a,w)
    op = exp(log(p) + ((a*z*sv)^2 - 2*a*v*z - (v^2)*x)/(2*(sv^2)*x+2))/sqrt((sv^2)*x+1)/(a^2);
end
end



function op = full_pdf(x,v,sv,a,z,sz,t,st,err)
%not sure if these are the actual defaults
n_st = 50;
n_sz = 50;
use_adaptive = 0;
simps_err = 1e-8;
%"""full pdf"""
%# Check if parpameters are valid (also don't compute if less than t)
if (z<0) || (z>1) || (a<0) || (t<0) || (st<0) || (sv<0) || (sz<0) || (sz>1) ||...
        ((abs(x)-(t-st/2.))<0) || (z+sz/2.>1) || (z-sz/2.<0) || (t-st/2.<0)
    op = 0;
    return;
end

%# transform x,v,z if x is upper bound response
if x > 0
    v = -v;
    z = 1.-z;
end
x = abs(x);

if st<1e-3
    st = 0;
end
if sz <1e-3
    sz = 0;
end

if (sz==0)
    if (st==0) %#sv=$,sz=0,st=0
        op = pdf_sv(x - t, v, sv, a, z, err);
    else      %#sv=$,sz=0,st=$
        if use_adaptive>0
            error('Not copied over');
        else
            op = simpson_1D(x, v, sv, a, z, t, err, z, z, 0, t-st/2., t+st/2., n_st);
        end
    end
else %#sz=$
    if (st==0) %#sv=0,sz=$,st=0
        if use_adaptive
            error('Not copied over');
        else
            op = simpson_1D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t, t , 0);
        end
    else      %#sv=0,sz=$,st=$
        if use_adaptive
            error('Not copied over');
        else
            op = simpson_2D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t-st/2., t+st/2., n_st);
        end
    end
end
end



function op = simpson_1D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
%assert ((n_sz&1)==0 and (n_st&1)==0), "n_st and n_sz have to be even"
%assert ((ub_t-lb_t)*(ub_z-lb_z)==0 and (n_sz*n_st)==0), "the function is defined for 1D-integration only"

n = max(n_st, n_sz);
if n_st==0 %integration over z
    hz = (ub_z-lb_z)/n;
    ht = 0;
    lb_t = t;
    ub_t = t;
else %integration over t
    hz = 0;
    ht = (ub_t-lb_t)/n;
    lb_z = z;
    ub_z = z;
end
S = pdf_sv(x - lb_t, v, sv, a, lb_z, err);


for i = 1:n
    z_tag = lb_z + hz * i;
    t_tag = lb_t + ht * i;
    y = pdf_sv(x - t_tag, v, sv, a, z_tag, err);
    if mod(i,2)==1 %check if i is odd
        S = S  + (4 * y);
    else
        S = S + (2 * y);
    end
end

S = S - y; %the last term should be f(b) and not 2*f(b) so we subtract y

S = S / ((ub_t-lb_t)+(ub_z-lb_z)); %the right function if pdf_sv()/sz or pdf_sv()/st

op = ((ht+hz) * S / 3);
end

function op = simpson_2D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
%assert ((n_sz&1)==0 and (n_st&1)==0), "n_st and n_sz have to be even"
%assert ((ub_t-lb_t)*(ub_z-lb_z)>0 and (n_sz*n_st)>0), "the function is defined for 2D-integration only, lb_t: %f, ub_t %f, lb_z %f, ub_z %f, n_sz: %d, n_st %d" % (lb_t, ub_t, lb_z, ub_z, n_sz, n_st)

%     cdef ht
%     cdef S
%     cdef t_tag, y
%     cdef i_t

ht = (ub_t-lb_t)/n_st;

S = simpson_1D(x, v, sv, a, z, lb_t, err, lb_z, ub_z, n_sz, 0, 0, 0);

for i_t  = 1:n_st
    t_tag = lb_t + ht * i_t;
    y = simpson_1D(x, v, sv, a, z, t_tag, err, lb_z, ub_z, n_sz, 0, 0, 0);
    if mod(i_t,2)==1 %check if i is odd
        S = S + (4 * y);
    else
        S = S + (2 * y);
    end
end
S = S - y; %the last term should be f(b) and not 2*f(b) so we subtract y
S = S/ (ub_t-lb_t);

op = (ht * S / 3);
end

