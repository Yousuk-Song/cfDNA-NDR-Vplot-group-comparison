#!/usr/bin/env python3
import argparse,csv,hashlib,itertools,json,subprocess,os
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.ndimage import gaussian_filter1d

def main():
 q=argparse.ArgumentParser();q.add_argument('--manifest',required=True);q.add_argument('--outdir',required=True);q.add_argument('--channel',default='STAD');q.add_argument('--sigma-bp',type=float,default=20);q.add_argument('--permutations',type=int,default=99999);q.add_argument('--seed',type=int,default=20261110);q.add_argument('--rscript',default='Rscript');a=q.parse_args();manifest=Path(a.manifest).resolve();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 with manifest.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 if not rows or not {'sample','group','path'}.issubset(rows[0]):raise SystemExit('TSV requires sample, group, path')
 if any(not r['sample'].strip() or not r['group'].strip() or not r['path'].strip() for r in rows) or len({r['sample'] for r in rows})!=len(rows):raise SystemExit('Empty field or duplicate sample')
 groups=list(dict.fromkeys(r['group'] for r in rows));
 if len(groups) not in (2,3):raise SystemExit('Exactly 2 or 3 groups required')
 if a.sigma_bp<0:raise SystemExit('sigma must be >=0')
 matrices=[];lengths=[];positions=[];audit=[];bedhash=set();unknown=[];edges=None
 for r in rows:
  p=Path(r['path']).expanduser();p=p if p.is_absolute() else manifest.parent/p;p=p/'vplot.npz' if p.is_dir() else p
  with np.load(p,allow_pickle=False) as z:
   raw=z['raw_counts'].astype(float);names=z['channels'].astype(str).tolist();xe=z['x_edges'];le=z['length_edges'];
   if raw.ndim!=3 or raw.shape!=(len(le)-1,len(xe)-1,len(names)) or not np.isfinite(raw).all() or (raw<0).any() or a.channel not in names:raise SystemExit(f'Invalid NPZ or channel: {p}')
   if edges is None:edges=(xe.copy(),le.copy())
   elif not(np.array_equal(xe,edges[0]) and np.array_equal(le,edges[1])):raise SystemExit('Bin edges differ')
   if not np.all(np.diff(xe)>0) or not np.allclose(np.diff(xe),np.diff(xe)[0]):raise SystemExit('Position bins must be uniform')
   lo,hi=raw.min(),raw.max();v=raw[:,:,names.index(a.channel)];total=v.sum()
   if hi<=lo or total<=0:raise SystemExit(f'Empty or constant matrix: {p}')
   matrices.append((v-lo)/(hi-lo));lengths.append(v.sum(axis=1)/total);positions.append(v.sum(axis=0)/total)
  meta=p.parent/'run.json';h=''
  if meta.exists():
   d=json.loads(meta.read_text());h=next((v for k,v in d.get('bed_sha256',{}).items() if k.split('=')[0]==a.channel),'')
  if h:bedhash.add(h)
  else:unknown.append(r['sample'])
  audit.append({'sample':r['sample'],'group':r['group'],'npz':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bed_sha256':h})
 if len(bedhash)>1:raise SystemExit('Different NDR BED hashes')
 M=np.stack(matrices);L=np.stack(lengths);P=np.stack(positions);
 if any(sum(r['group']==g for r in rows)<2 for g in groups):raise SystemExit('At least two patients per group are required')
 if 'central_coverage' not in rows[0] and (out/'central_coverage_test.tsv').exists(): (out/'central_coverage_test.tsv').unlink()
 np.savetxt(out/'patient_length_profiles.tsv',L,delimiter='\t');np.savetxt(out/'patient_midpoint_profiles.tsv',P,delimiter='\t'); pd.DataFrame(rows).to_csv(out/'patient_groups.tsv',sep='\t',index=False)
 subprocess.run([a.rscript,str(Path(__file__).with_name('statistics.R')),str(out/'patient_length_profiles.tsv'),str(out/'patient_midpoint_profiles.tsv'),str(out/'patient_groups.tsv'),str(out),str(a.permutations),str(a.seed)],check=True)
 tests=pd.read_csv(out/'distribution_tests.tsv',sep='\t').set_index('feature');central_path=out/'central_coverage_test.tsv'; central=pd.read_csv(central_path,sep='\t').iloc[0] if 'central_coverage' in rows[0] else None;idx=[np.array([r['group']==g for r in rows]) for g in groups];means=[M[i].mean(axis=0) for i in idx];pairs=list(itertools.combinations(range(len(groups)),2));diffs=[means[j]-means[i] for i,j in pairs];lim=max(float(np.abs(d).max()) for d in diffs) or 1.;cols=['#0072B2','#009E73','#D55E00'];nplot=len(groups)+len(pairs);fig=plt.figure(figsize=(4*nplot+10,6.8),layout='constrained');gs=fig.add_gridspec(1,nplot+2,width_ratios=[1]*nplot+[1,1.35]);original_axes=[fig.add_subplot(gs[0,j]) for j in range(len(groups))];axes=[fig.add_subplot(gs[0,len(groups)+j]) for j in range(len(pairs))];al=fig.add_subplot(gs[0,-2]);ap=fig.add_subplot(gs[0,-1]);xe,le=edges;x=(xe[:-1]+xe[1:])/2;length=(le[:-1]+le[1:])/2
 vmax=max(float(v.max()) for v in means) or 1.
 for ax,g,v,ii in zip(original_axes,groups,means,idx):
  original_im=ax.imshow(v,origin='lower',extent=[xe[0],xe[-1],le[0],le[-1]],aspect='auto',cmap='viridis',vmin=0,vmax=vmax,interpolation='nearest');ax.axvline(0,color='white',ls=':',lw=.7);ax.set(title=f'{g} (n={ii.sum()})',xlabel='Distance from NDR midpoint (bp)',ylabel='Fragment length (bp)')
 fig.colorbar(original_im,ax=original_axes,label='Mean min–max count',shrink=.8)
 for ax,(i,j),d in zip(axes,pairs,diffs):
  im=ax.imshow(d,origin='lower',extent=[xe[0],xe[-1],le[0],le[-1]],aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),interpolation='nearest');ax.axvline(0,color='black',ls=':',lw=.7);ax.set(title=f'{groups[j]} − {groups[i]}',xlabel='Distance from NDR midpoint (bp)',ylabel='Fragment length (bp)')
 fig.colorbar(im,ax=axes,label='Difference in mean min–max count',shrink=.8)
 for g,col,ii in zip(groups,cols,idx):
  v=L[ii].mean(0);w=P[ii].mean(0);al.plot(v,length,color=col,lw=2,label=f'{g} (n={ii.sum()})');al.fill_betweenx(length,0,v,color=col,alpha=.12);sigma=a.sigma_bp/float(np.diff(xe)[0]);ap.plot(x,gaussian_filter1d(w,sigma,mode='reflect') if sigma else w,color=col,lw=2,label=f'{g} (n={ii.sum()})')
 al.set(title='B. Fragment-length distribution',xlabel='Proportion per length bin',ylabel='Fragment length (bp)',xlim=(0,None),ylim=(le[0],le[-1]));ap.set(title='C. Fragment-midpoint distribution',xlabel='Distance from NDR midpoint (bp)',ylabel='Proportion per position bin',xlim=(xe[0],xe[-1]));al.legend(frameon=False,loc='lower right',fontsize=9);ap.axvline(0,color='gray',ls=':',lw=.7);ap.legend(frameon=False,loc='lower left');original_axes[0].text(0,1.165,'A. Group V-plots and differences',transform=original_axes[0].transAxes,ha='left',fontsize=12,fontweight='bold');al.set_title('B. Fragment-length distribution',y=1.15,fontsize=12,fontweight='bold');ap.set_title('C. Fragment-midpoint distribution',y=1.15,fontsize=12,fontweight='bold');al.text(.5,1.025,'PERMANOVA P = %.4g (q = %.4g) | PERMDISP P = %.4g' % (tests.loc['fragment_length','P'],tests.loc['fragment_length','BH_q_BC'],tests.loc['fragment_length','PERMDISP_P']),transform=al.transAxes,ha='center',va='bottom',fontsize=9,fontweight='normal');ap.text(.5,1.025,'PERMANOVA P = %.4g (q = %.4g) | PERMDISP P = %.4g' % (tests.loc['midpoint','P'],tests.loc['midpoint','BH_q_BC'],tests.loc['midpoint','PERMDISP_P']),transform=ap.transAxes,ha='center',va='bottom',fontsize=9,fontweight='normal');fig.suptitle('NDR V-plot group comparison' + (f" | Griffin central ±30 bp permutation P = {central['P']:.4g}" if central is not None else ''),fontsize=15);fig.savefig(out/'Vplot_group_comparison.png',dpi=230,bbox_inches='tight');fig.savefig(out/'Vplot_group_comparison.pdf',bbox_inches='tight');plt.close(fig)
 with (out/'input_audit.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=audit[0],delimiter='\t');w.writeheader();w.writerows(audit)
 np.savez_compressed(out/'group_matrices.npz',groups=np.array(groups),means=np.stack(means),differences=np.stack(diffs),difference_names=np.array([f'{groups[j]} - {groups[i]}' for i,j in pairs]),length_proportions=np.stack([L[ii].mean(0) for ii in idx]),midpoint_proportions=np.stack([P[ii].mean(0) for ii in idx]),x_edges=xe,length_edges=le)
 (out/'analysis_metadata.json').write_text(json.dumps({'group_counts':{g:int(ii.sum()) for g,ii in zip(groups,idx)},'channel':a.channel,'difference_scale':lim,'sigma_bp':a.sigma_bp,'bed_hashes':sorted(bedhash),'bed_hash_unverified_samples':unknown,'normalization_A':'per-sample minmax across all channels','normalization_BC':'per-sample channel raw counts / total','statistical_test':'Euclidean PERMANOVA and spatial-median PERMDISP, unsmoothed patient proportions; BH over B and C','permutations':a.permutations,'seed':a.seed},indent=2)+'\n')
 if unknown:print('WARNING: BED provenance unverified for',len(unknown),'samples')
 print('DONE',out)
if __name__=='__main__':main()
