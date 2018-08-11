function op = hddm_pdf_full(x,v,sv,a,z,sz,t,st,err)

%Translated from HDDM core
n_st = 15;
n_sz = 15;
use_adaptive = 0;
%             simps_err = 1e-8;
%"""full pdf"""
%% Check if parpameters are valid (also don't compute if less than t)
if (z<0) || (z>1) || (a<0) || (t<0) || (st<0) || (sv<0) || (sz<0) || (sz>1) ||...
        ((abs(x)-(t-st/2.))<0) || (z+sz/2.>1) || (z-sz/2.<0) || (t-st/2.<0)
    op = 0;
    return;
end

%% transform x,v,z if x is upper bound response
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
    if (st==0) %%sv=$,sz=0,st=0
        op = hddm_pdf_sv(x - t, v, sv, a, z, err);
    else      %%sv=$,sz=0,st=$
        if use_adaptive>0
            error('Not copied over');
        else
            op = hddm_simpson_1D(x, v, sv, a, z, t, err, z, z, 0, t-st/2., t+st/2., n_st);
        end
    end
else %%sz=$
    if (st==0) %%sv=0,sz=$,st=0
        if use_adaptive
            error('Not copied over');
        else
            op = hddm_simpson_1D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t, t , 0);
        end
    else      %%sv=0,sz=$,st=$
        if use_adaptive
            error('Not copied over');
        else
            op = hddm_simpson_2D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t-st/2., t+st/2., n_st);
        end
    end
end
end

function op = hddm_pdf(x,v,a,w,err)
%Translated from HDDM core - this is so that fast fitting can
%be done when we don't need to use extended models. Or when we
%want to find some good initisalisation parameters for extended
%models.
%"""Compute the likelihood of the drift diffusion model f(t|v,a,z) using the method
%and implementation of Navarro & Fuss, 2009.
%"""
if x <= 0
    op = 0;
else
    
    tt = x/(a^2);% # use normalized time
    p = navfuss_ftt_01w(tt, w, err);% #get f(t|0,1,w)
    
    %% convert to f(t|v,a,w)
    op =  p*exp(-v*a*w -(v^2)*x/2)/(a^2);
end
end

function op = hddm_pdf_sv(x,v,sv,a,z,err)
%Translated from HDDM core
%"""Compute the likelihood of the drift diffusion model f(t|v,a,z,sv) using the method
%and implementation of Navarro & Fuss, 2009.
%sv is the std of the drift rate
%"""
if x <= 0
    op = 0;
elseif sv==0
    op =  hddm_pdf(x, v, a, z, err);
else
    tt = x/(a^2); %% use normalized time
    p  = navfuss_ftt_01w(tt, z, err); %%get f(t|0,1,w)
    
    %% convert to f(t|v,a,w)
    op = exp(log(p) + ((a*z*sv)^2 - 2*a*v*z - (v^2)*x)/(2*(sv^2)*x+2))/sqrt((sv^2)*x+1)/(a^2);
end
end

function op = hddm_simpson_1D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
%Translated from HDDM core
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
S = hddm_pdf_sv(x - lb_t, v, sv, a, lb_z, err);


for i = 1:n
    z_tag = lb_z + hz * i;
    t_tag = lb_t + ht * i;
    y = hddm_pdf_sv(x - t_tag, v, sv, a, z_tag, err);
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



function op = hddm_simpson_2D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
%Translated from HDDM core
%assert ((n_sz&1)==0 and (n_st&1)==0), "n_st and n_sz have to be even"
%assert ((ub_t-lb_t)*(ub_z-lb_z)>0 and (n_sz*n_st)>0), "the function is defined for 2D-integration only, lb_t: %f, ub_t %f, lb_z %f, ub_z %f, n_sz: %d, n_st %d" % (lb_t, ub_t, lb_z, ub_z, n_sz, n_st)

ht = (ub_t-lb_t)/n_st;
S = hddm_simpson_1D(x, v, sv, a, z, lb_t, err, lb_z, ub_z, n_sz, 0, 0, 0);

for i_t  = 1:n_st
    t_tag = lb_t + ht * i_t;
    y = hddm_simpson_1D(x, v, sv, a, z, t_tag, err, lb_z, ub_z, n_sz, 0, 0, 0);
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