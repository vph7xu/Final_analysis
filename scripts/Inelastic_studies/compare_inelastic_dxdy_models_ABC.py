#!/usr/bin/env python3
"""GEn-II GEN3 helicity-summed 2D template fits:

    M(dx,dy) = A * [ P(dx,dy) + B*N(dx,dy) + C*I(dx,dy) ]

P/N: QE SIMC unit-normalized 2D templates (fnucl=1/0).
I: simulated inelastic template OR nonnegative tensor-product Bernstein
   polynomial of degree n in *each* of dx and dy (Pol2..Pol5).
All shapes have unit bin sum in the 2D fitting domain.
A = proton yield; B = neutron/proton; C = inelastic/proton.

Fitted with binned Poisson deviance, *all* cells including empty cells.
Reports Pearson reduced chi2, deviance/ndf and AIC; the chi2 p-value
is only approximate for very sparse two-dimensional histograms.

QE and simulated inelastic templates are accumulated on a WIDER 2D grid
than the fit/data grid. Shifts use integrated 2D histogram bins (a bilinear
CDF), then each component is normalized over the UNCHANGED fit region.
This avoids artificial template discontinuities at the fit boundaries.
NOTE: no manual GEN3 pre-shifts by default; set them explicitly if desired.
See --help for the switches that float/fix 2D template offsets.
"""
from __future__ import annotations
import argparse
import csv
import json
from math import comb
from pathlib import Path
import numpy as np

DATA = '/media/vimukthi/Sata/Research_JLab/outfiles/pass3/intimehighest/QE_data_GEN3_sbs100p_nucleon_np_model2_sbstrackingon_intimehighest.root'
QE = '/media/vimukthi/Sata/Research_JLab/outfiles/SIMC/noNRad_correct_raster/QE_sim_simc_He3eN_GEN3_sbs100p_nucleon_np_model2_elastic_simc.root'
INEL = '/media/vimukthi/Sata/Research_JLab/outfiles/inelastic_simulation_2025/QE_sim_GEN3_sbs100p_nucleon_np_model2_inelastic.root'


def arguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--data', default=DATA)
    p.add_argument('--qe', default=QE)
    p.add_argument('--inel', default=INEL)
    p.add_argument('--cuts', default='cuts/GEN3.json')
    p.add_argument('--run-quality', default='DB/run_quality.json')
    p.add_argument('--ignore-run-quality', action='store_true')
    p.add_argument('--data-tree', default=None)
    p.add_argument('--qe-tree', default=None)
    p.add_argument('--inel-tree', default=None)
    p.add_argument('--bins-dx', type=int, default=60)
    p.add_argument('--bins-dy', type=int, default=40)
    p.add_argument('--dx-range', nargs=2, type=float, default=(-2.2, 0.8))
    p.add_argument('--dy-padding', type=float, default=0.05)
    p.add_argument('--degrees', nargs='+', type=int, default=(2,3,4,5))
    p.add_argument('--float-shifts', action='store_true', help='Fit dx/dy shifts for P,N and simulated I; polynomial I has no shift')
    p.add_argument('--shift-x-max', type=float, default=0.05)
    p.add_argument('--shift-y-max', type=float, default=0.03)
    p.add_argument('--template-guard-dx', type=float, default=0.45, help='Minimum extra MC support on EACH dx side [m]')
    p.add_argument('--template-guard-dy', type=float, default=0.15, help='Minimum extra MC support on EACH dy side [m]')
    p.add_argument('--pre-shift-p-dx', type=float, default=0., help='Manual pre-shift in dx of QE proton (GEN3 original 2D: +0.05)')
    p.add_argument('--pre-shift-n-dx', type=float, default=0.)
    p.add_argument('--pre-shift-i-dx', type=float, default=0., help='Manual pre-shift in dx of inelastic (GEN3 original 2D: -0.10)')
    p.add_argument('--max-events', type=int, default=None, help='Max entries PER ROOT file (debug)')
    p.add_argument('--step-size', default='100 MB')
    p.add_argument('--maxiter', type=int, default=1000)
    p.add_argument('--outdir', default='out_plots/inelastic_dxdy_background_models_GEN3_ABC_py')
    return p.parse_args()


def good_runs(filename):
    with open(filename, encoding='utf-8') as f:
        j = json.load(f)
    runs = set()
    for key, val in j.items():
        if not isinstance(val, dict):
            continue
        if val.get('helicity') in (True, 1) and val.get('moller') in (True, 1):
            try: runs.add(int(key))
            except ValueError: pass
    if not runs:
        raise ValueError(f'No good runs found in {filename}; inspect run-quality JSON or use --ignore-run-quality')
    return runs


def tree_name(file, override):
    if override:
        if override not in file: raise KeyError(f'Missing tree {override} in {file.file_path}')
        return override
    for name in ('Tout','T'):
        if name in file: return name
    raise KeyError(f'Could not locate Tout or T in {file.file_path}')


def fill_histograms(args, cuts, runs):
    import uproot
    nx, ny = args.bins_dx, args.bins_dy
    xedges = np.linspace(*args.dx_range, nx + 1)
    yedges = np.linspace(cuts['dy_L']-args.dy_padding, cuts['dy_H']+args.dy_padding, ny+1)
    wx, wy = xedges[1]-xedges[0], yedges[1]-yedges[0]

    # Require enough extra support for ALL allowed floating shifts and at
    # least three additional bins. Use integer guard bins to align axes.
    gx = int(np.ceil(max(args.template_guard_dx, args.shift_x_max + 3*wx) / wx))
    gy = int(np.ceil(max(args.template_guard_dy, args.shift_y_max + 3*wy) / wy))
    source_xedges = xedges[0] + np.arange(-gx, nx+gx+1, dtype=float)*wx
    source_yedges = yedges[0] + np.arange(-gy, ny+gy+1, dtype=float)*wy
    print(f'[support] fit dx=[{xedges[0]:.4f},{xedges[-1]:.4f}] dy=[{yedges[0]:.4f},{yedges[-1]:.4f}]',flush=True)
    print(f'[support] MC  dx=[{source_xedges[0]:.4f},{source_xedges[-1]:.4f}] '
          f'dy=[{source_yedges[0]:.4f},{source_yedges[-1]:.4f}] '
          f'(guard {gx} x-bins, {gy} y-bins)',flush=True)

    h = {'data': np.zeros((nx,ny),dtype='float64')}
    mc = {k: np.zeros((nx+2*gx, ny+2*gy),dtype='float64')
          for k in ('proton','neutron','inelastic')}
    common = ['dx','dy','vz','eHCAL','ePS','eSH','trP','W2']
    inputs = [('data',args.data,args.data_tree),('qe',args.qe,args.qe_tree),('inel',args.inel,args.inel_tree)]
    for kind, filename, override in inputs:
        extra = ['runnum','ntrack','coin_time','helicity'] if kind=='data' else (['weight','fnucl'] if kind=='qe' else ['weight'])
        needed = common+extra
        nread = nselected = 0
        with uproot.open(filename) as root:
            name = tree_name(root,override)
            tree = root[name]
            missing = set(needed) - set(tree.keys())
            if missing: raise KeyError(f'{filename}:{name} missing branches: {sorted(missing)}')
            print(f'[input] {kind:>6s}: {filename} ({tree.num_entries:,} entries)',flush=True)
            for block in tree.iterate(expressions=needed,library='np',step_size=args.step_size):
                if not isinstance(block,dict):
                    block={key:np.asarray(block[key]) for key in needed}
                count=len(block['dx'])
                if args.max_events is not None:
                    room=args.max_events-nread
                    if room<=0: break
                    if count>room:
                        block={k:v[:room] for k,v in block.items()}
                        count=room
                nread+=count
                with np.errstate(divide='ignore',invalid='ignore'):
                    ep=(block['ePS']+block['eSH'])/block['trP']
                # Events must pass the SAME non-geometric physics cuts; the
                # wider MC dx/dy support is only an interpolation reservoir.
                dylo,dyhi=(yedges[0],yedges[-1]) if kind=='data' else (source_yedges[0],source_yedges[-1])
                ok=(np.isfinite(block['dx']) & np.isfinite(block['dy']) & np.isfinite(ep) &
                    (np.abs(block['vz']) < 0.27) & (block['eHCAL']>cuts['eHCAL_L']) &
                    (block['ePS']>0.2) & (ep>cuts['EoverP_L']) & (ep<cuts['EoverP_H']) &
                    (block['W2']>cuts['W2_L']) & (block['W2']<cuts['W2_H']) &
                    (block['dy']>dylo) & (block['dy']<dyhi))
                if kind=='data':
                    ok &= ((block['runnum']>=cuts['runnum_L']) & (block['runnum']<=cuts['runnum_H']) &
                           (block['ntrack']>=1) & (np.abs(block['helicity'])==1) &
                           (block['coin_time']>=cuts['coin_L']) & (block['coin_time']<=cuts['coin_H']))
                    if runs is not None:
                        ok &= np.isin(block['runnum'],list(runs))
                    vals=np.column_stack((block['dx'][ok],block['dy'][ok]))
                    h['data']+=np.histogramdd(vals,bins=[xedges,yedges])[0]
                elif kind=='qe':
                    ok &= np.isfinite(block['weight'])
                    for nucl,key,preshift in ((1,'proton',args.pre_shift_p_dx),(0,'neutron',args.pre_shift_n_dx)):
                        sel=ok & (block['fnucl']==nucl)
                        vals=np.column_stack((block['dx'][sel]+preshift,block['dy'][sel]))
                        mc[key]+=np.histogramdd(vals,bins=[source_xedges,source_yedges],weights=block['weight'][sel])[0]
                else:
                    ok &= np.isfinite(block['weight'])
                    vals=np.column_stack((block['dx'][ok]+args.pre_shift_i_dx,block['dy'][ok]))
                    mc['inelastic']+=np.histogramdd(vals,bins=[source_xedges,source_yedges],weights=block['weight'][ok])[0]
                nselected+=int(np.count_nonzero(ok))
                if args.max_events is not None and nread>=args.max_events: break
        print(f'        passed non-dx cuts: {nselected:,} / {nread:,}',flush=True)
    for name,v in {**h,**mc}.items():
        if v.sum()<=0 or not np.isfinite(v).all() or (v<-1e-9*abs(v.sum())).any():
            raise ValueError(f'Invalid histogram {name} (sum={v.sum()}); check cuts and MC weight signs')
    return h,mc,xedges,yedges,source_xedges,source_yedges


class ShiftedTemplate2D:
    """Bin-integrated translated MC histogram using its bilinear 2D CDF.

    The MC is histogrammed on extended support. Each fit bin receives the
    integral of the piecewise-constant original density over the rectangle
    (xlo-sx,xhi-sx) x (ylo-sy,yhi-sy). The result is normalized in FIT bins.
    """
    def __init__(self, array, source_xedges, source_yedges, fit_xedges, fit_yedges):
        self.sx=np.asarray(source_xedges)
        self.sy=np.asarray(source_yedges)
        self.fx=np.asarray(fit_xedges)
        self.fy=np.asarray(fit_yedges)
        self.wx=self.sx[1]-self.sx[0]
        self.wy=self.sy[1]-self.sy[0]
        a=np.maximum(np.asarray(array,dtype=float),0.0)
        if a.sum()<=0: raise ValueError('Empty source MC template')
        a=a/a.sum()
        self.cdf=np.pad(a.cumsum(axis=0).cumsum(axis=1),((1,0),(1,0)))
        self.nx,self.ny=a.shape

    def shifted(self, sx=0., sy=0.):
        # Bilinear interpolation in CDF coordinates; clipping corresponds to
        # zero density outside the extended MC support.
        u=np.clip((self.fx-sx-self.sx[0])/self.wx,0.,self.nx)
        v=np.clip((self.fy-sy-self.sy[0])/self.wy,0.,self.ny)
        i=np.minimum(np.floor(u).astype(int),self.nx-1)
        j=np.minimum(np.floor(v).astype(int),self.ny-1)
        a=u-i; b=v-j
        c=self.cdf
        f=(c[i[:,None],j[None,:]]*(1-a[:,None])*(1-b[None,:])+
           c[(i+1)[:,None],j[None,:]]*a[:,None]*(1-b[None,:])+
           c[i[:,None],(j+1)[None,:]]*(1-a[:,None])*b[None,:]+
           c[(i+1)[:,None],(j+1)[None,:]]*a[:,None]*b[None,:])
        masses=np.maximum(np.diff(np.diff(f,axis=0),axis=1),0.)
        total=masses.sum()
        if total<=0: raise ValueError('Shifted MC has no support within 2D fit domain')
        return masses/total


def bernstein_1d(n, centers, low, high):
    t=np.clip((centers-low)/(high-low),0,1)
    cols=np.stack([comb(n,k)*t**k*(1-t)**(n-k) for k in range(n+1)],axis=1)
    return cols/cols.sum(axis=0)[None,:]


def poly_shape(degree, pars, basis_x, basis_y):
    """Positive 2D polynomial with unit bin sum; degree n in each axis."""
    k=(degree+1)**2
    z=np.r_[pars,np.zeros(1)]
    z=np.exp(z-z.max())
    weights=(z/z.sum()).reshape(degree+1,degree+1)
    out=basis_x@weights@basis_y.T
    return out/out.sum()


class FitProblem:
    def __init__(self,h,mc,xedges,yedges,source_xedges,source_yedges,degree,float_shifts):
        self.data=h['data']
        self.templates={k:ShiftedTemplate2D(v,source_xedges,source_yedges,xedges,yedges)
                        for k,v in mc.items()}
        self.degree=degree
        self.float_shifts=float_shifts
        self.xedges=xedges; self.yedges=yedges
        self.wx=xedges[1]-xedges[0]; self.wy=yedges[1]-yedges[0]
        self.nshape=0 if degree<0 else (degree+1)**2-1
        self.nshifts=0 if not float_shifts else (6 if degree<0 else 4)
        self.npar=3+self.nshape+self.nshifts
        self.nobs=self.data.sum()
        if degree>=0:
            self.bx=bernstein_1d(degree,0.5*(xedges[:-1]+xedges[1:]),xedges[0],xedges[-1])
            self.by=bernstein_1d(degree,0.5*(yedges[:-1]+yedges[1:]),yedges[0],yedges[-1])

    def pieces(self,par):
        a,b,c=par[:3]
        ix=3+self.nshape
        if self.float_shifts:
            sxp,syp,sxn,syn=par[ix:ix+4]
            sxi,syi=(par[ix+4:ix+6] if self.degree<0 else (0,0))
        else: sxp=syp=sxn=syn=sxi=syi=0.
        p=self.templates['proton'].shifted(sxp,syp)
        n=self.templates['neutron'].shifted(sxn,syn)
        if self.degree<0:
            I=self.templates['inelastic'].shifted(sxi,syi)
        else:
            I=poly_shape(self.degree,par[3:3+self.nshape],self.bx,self.by)
        return a*p,a*b*n,a*c*I

    def predict(self,par):
        p,n,I=self.pieces(par)
        return p+n+I

    def deviance(self,par):
        mu=np.maximum(self.predict(par),1e-10)
        n=self.data
        terms=mu-n
        select=n>0
        terms[select]+=n[select]*np.log(n[select]/mu[select])
        return float(2*terms.sum())


def fit_model(problem,args,label):
    from scipy.optimize import minimize
    from scipy.stats import chi2 as chi2dist
    n=problem.nobs
    nsh=problem.nshape
    # Work in dimensionless proton normalization a=A/Ndata. Without this
    # rescaling, finite-difference L-BFGS-B often leaves A close to its
    # INITIAL value while incorrectly compensating with B and C.
    start=np.r_[1/1.5,0.2,0.3,np.zeros(nsh),np.zeros(problem.nshifts)]
    bounds=[(1e-10,3.0),(0,20),(0,20)]+[(-9,9)]*nsh
    if problem.float_shifts:
        bounds+=([(-args.shift_x_max,args.shift_x_max),(-args.shift_y_max,args.shift_y_max)]*
                 (3 if problem.degree<0 else 2))
    def physical_parameters(theta):
        p=theta.copy()
        p[0]*=n
        return p
    # Large high-order polynomial models are strongly correlated; L-BFGS-B first.
    result=minimize(lambda theta:problem.deviance(physical_parameters(theta)),
                    start,method='L-BFGS-B',bounds=bounds,
                    options={'maxiter':args.maxiter,'ftol':5e-10,'maxls':40})
    par=physical_parameters(result.x)
    pieces=problem.pieces(par)
    mu=np.maximum(sum(pieces),1e-10)
    obs=problem.data
    pearson=np.sum((obs-mu)**2/mu)
    dev=problem.deviance(par)
    dof=obs.size-problem.npar
    if dof<=0: raise ValueError('More parameters than bins: increase --bins-dx/--bins-dy')
    x=0.5*(problem.xedges[:-1]+problem.xedges[1:])
    y=0.5*(problem.yedges[:-1]+problem.yedges[1:])
    return dict(model=label,success=result.success,status=str(result.message),par=par,
                pieces=pieces,fit=mu,resid=obs-mu,pull=(obs-mu)/np.sqrt(mu),
                chi2=float(pearson),ndf=int(dof),reduced=float(pearson/dof),
                deviance=dev,deviance_reduced=dev/dof,p_value=float(chi2dist.sf(pearson,dof)),
                aic=dev+2*problem.npar,x=x,y=y,n_parameters=problem.npar)


def plot_fit(r,h,edges,outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Ellipse
    x,y=r['x'],r['y']
    xx,yy=edges
    proton,neutron,inel=r['pieces']
    arrays=[h['data'],r['fit'],r['resid'],r['pull']]
    titles=['Data','Total fit','Data - fit','Pearson pull']
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for a,ar,title in zip(ax.flat,arrays,titles):
        kwargs={}
        if title in ('Data - fit','Pearson pull'):
            largest=np.quantile(np.abs(ar),0.99)
            kwargs['vmin']=-max(largest,1e-9);kwargs['vmax']=max(largest,1e-9)
            kwargs['cmap']='RdBu_r'
        im=a.pcolormesh(xx,yy,ar.T,shading='auto',**kwargs)
        a.set_xlabel('dx [m]');a.set_ylabel('dy [m]');a.set_title(title)
        fig.colorbar(im,ax=a)
    fig.suptitle(f"{r['model']}: Pearson $\\chi^2/\\mathrm{{ndf}}$ = {r['chi2']:.1f}/{r['ndf']} = {r['reduced']:.3f}")
    fig.savefig(outdir/f"fit2d_{r['model'].lower()}.pdf")
    fig.savefig(outdir/f"fit2d_{r['model'].lower()}.png",dpi=160)
    plt.close(fig)
    fig,(a,b)=plt.subplots(2,2,figsize=(12,8),layout='constrained',gridspec_kw={'height_ratios':[3,1]})
    widths=[np.diff(xx)[0],np.diff(yy)[0]]
    for col,(centers,axislabel,axis) in enumerate(((x,'dx [m]',1),(y,'dy [m]',0))):
        top,bot=a[col],b[col]
        data=h['data'].sum(axis=axis)
        p=proton.sum(axis=axis);n=neutron.sum(axis=axis);ine=inel.sum(axis=axis)
        model=r['fit'].sum(axis=axis)
        top.errorbar(centers,data,yerr=np.sqrt(np.maximum(data,0)),fmt='.',markersize=4,color='k',label='Data')
        for z,label in ((model,'Total'),(p,'QE p'),(n,'QE n'),(ine,'Background')):
            top.step(centers,z,where='mid',label=label)
        top.set_ylabel('Counts / bin');top.legend(fontsize=8,ncol=2)
        residual=data-model
        pull=residual/np.sqrt(np.maximum(model,1e-10))
        bot.axhline(0,color='0.3',linestyle='--',lw=1)
        bot.plot(centers,pull,'.',markersize=4)
        bot.set_xlabel(axislabel);bot.set_ylabel('Pull')
        bound=max(4,np.nanquantile(np.abs(pull),0.98)*1.2)
        bot.set_ylim(-bound,bound)
    fig.suptitle(f"{r['model']}  |  A={r['par'][0]:.1f}, B={r['par'][1]:.3f}, C={r['par'][2]:.3f}")
    fig.savefig(outdir/f"projections_{r['model'].lower()}.pdf")
    fig.savefig(outdir/f"projections_{r['model'].lower()}.png",dpi=160)
    plt.close(fig)


def run(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if args.bins_dx < 8 or args.bins_dy < 8: raise ValueError('Need at least eight bins along each axis')
    if max(args.degrees) > 8 or min(args.degrees) < 0: raise ValueError('Use polynomial degrees 0..8')
    if args.dx_range[0]>=args.dx_range[1]: raise ValueError('Invalid dx range')
    if args.template_guard_dx<0 or args.template_guard_dy<0: raise ValueError('Template guards must be nonnegative')
    if args.shift_x_max<0 or args.shift_y_max<0: raise ValueError('Shift bounds must be nonnegative')
    with open(args.cuts,encoding='utf-8') as f: cuts=json.load(f)
    for key in ('dy_L','dy_H','dx_L','dx_H','dx_r','dy_r'):
        if key not in cuts: raise KeyError(f'Missing {key} in {args.cuts}')
    runs=None if args.ignore_run_quality else good_runs(args.run_quality)
    h,mc,xe,ye,source_xe,source_ye=fill_histograms(args,cuts,runs)
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True)
    x=0.5*(xe[1:]+xe[:-1]);y=0.5*(ye[1:]+ye[:-1])
    DX,DY=np.meshgrid(x,y,indexing='ij')
    ellipse=(DX/cuts['dx_r'])**2+(DY/cuts['dy_r'])**2<=1.
    dx_window=(DX>=cuts['dx_L'])&(DX<=cuts['dx_H'])
    fits=[]
    for name,degree in [('Simulation',-1)]+[(f'Pol{k}',k) for k in sorted(set(args.degrees))]:
        pr=FitProblem(h,mc,xe,ye,source_xe,source_ye,degree,args.float_shifts)
        print(f'[fit] {name}: {pr.npar} parameters, {h["data"].size} bins',flush=True)
        r=fit_model(pr,args,name)
        r['bg_frac_ellipse']=float(r['pieces'][2][ellipse].sum()/h['data'][ellipse].sum()) if h['data'][ellipse].sum()>0 else np.nan
        r['bg_frac_dx_window']=float(r['pieces'][2][dx_window].sum()/h['data'][dx_window].sum()) if h['data'][dx_window].sum()>0 else np.nan
        print(f"      converged={r['success']} chi2/ndf={r['reduced']:.3f}  deviance/ndf={r['deviance_reduced']:.3f} A={r['par'][0]:.3g} B={r['par'][1]:.3g} C={r['par'][2]:.3g} f_bg(ellipse)={r['bg_frac_ellipse']:.4f}",flush=True)
        if not r['success']: print(f"[warning] optimizer: {r['status']}")
        fits.append(r)
        plot_fit(r,h,(xe,ye),outdir)
    columns=['model','success','A','B','C','n_parameters','ndf','pearson_chi2','reduced_chi2','deviance','deviance_per_ndf','p_value_approx','AIC','inelastic_counts_full','inelastic_counts_ellipse','data_counts_ellipse','bg_frac_ellipse','bg_frac_dx_window']
    with open(outdir/'fit_summary.csv','w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for r in fits:
            bg=r['pieces'][2]
            writer.writerow(dict(model=r['model'],success=r['success'],A=r['par'][0],B=r['par'][1],C=r['par'][2],
                                 n_parameters=r['n_parameters'],ndf=r['ndf'],pearson_chi2=r['chi2'],reduced_chi2=r['reduced'],
                                 deviance=r['deviance'],deviance_per_ndf=r['deviance_reduced'],p_value_approx=r['p_value'],AIC=r['aic'],
                                 inelastic_counts_full=bg.sum(),inelastic_counts_ellipse=bg[ellipse].sum(),data_counts_ellipse=h['data'][ellipse].sum(),
                                 bg_frac_ellipse=r['bg_frac_ellipse'],bg_frac_dx_window=r['bg_frac_dx_window']))
    with open(outdir/'fit_parameters.csv','w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['model','parameter','value'])
        for r in fits:
            for k,v in enumerate(r['par']):
                pname=['A','B','C'][k] if k<3 else f'shape_or_shift_param_{k}'
                writer.writerow([r['model'],pname,v])
    np.savez_compressed(outdir/'fit_arrays.npz',dx_edges=xe,dy_edges=ye,data=h['data'],
                        **{f'{r["model"]}_{key}':value for r in fits for key,value in
                           [('model',r['fit']),('p',r['pieces'][0]),('n',r['pieces'][1]),('inel',r['pieces'][2]),('residual',r['resid']),('pull',r['pull'])]})
    labels=[r['model'] for r in fits]
    positions=np.arange(len(fits))
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    ax.bar(positions,[r['reduced'] for r in fits]);ax.set_xticks(positions,labels)
    ax.axhline(1,color='k',ls='--',lw=1);ax.set_ylabel('Pearson $\\chi^2$/ndf')
    ax.set_title('2D model comparison (all selected bins)')
    fig.savefig(outdir/'reduced_chi2_comparison.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    ax.plot(positions,[r['bg_frac_ellipse'] for r in fits],'o-')
    ax.set_xticks(positions,labels);ax.set_ylabel('Fitted inelastic / data in neutron ellipse')
    ax.set_title('2D neutron-ellipse background fraction (bin-center selection)')
    fig.savefig(outdir/'background_fraction_ellipse.pdf');plt.close(fig)
    with open(outdir/'configuration.json','w') as f:
        json.dump({**vars(args),'model':'A*(P+B*N+C*I)','poly':'tensor-product Bernstein, degree n in each coordinate; unit-bin-sum',
                   'ellipse':'(dx/dx_r)^2+(dy/dy_r)^2 <= 1; bin-center inclusion',
                   'source_dx_edges':[float(source_xe[0]),float(source_xe[-1])],
                   'source_dy_edges':[float(source_ye[0]),float(source_ye[-1])],
                   'shifted_template':'2D bin-integrated CDF on extended MC support; normalize inside fit grid',
                   'cuts_roi':{k:cuts[k] for k in ('dx_L','dx_H','dx_r','dy_r','dy_L','dy_H')}},f,indent=2)
    print(f'[done] {outdir}')

if __name__=='__main__':
    run(arguments())
