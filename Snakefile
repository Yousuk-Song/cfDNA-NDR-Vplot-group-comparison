from pathlib import Path
import csv
BASE=Path(workflow.basedir)
manifest=Path(config['manifest']).expanduser().resolve()
out=Path(config['output_dir']).expanduser().resolve()
with manifest.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
paths=[]
for row in rows:
 p=Path(row['path']).expanduser();p=p if p.is_absolute() else manifest.parent/p
 paths.append(str((p/'vplot.npz' if p.is_dir() else p).resolve()))
output_names=(['central_coverage_test.tsv'] if rows and 'central_coverage' in rows[0] else []) + ['Vplot_group_comparison.png','Vplot_group_comparison.pdf','input_audit.tsv','group_matrices.npz','analysis_metadata.json','distribution_tests.tsv','R_sessionInfo.txt']
rule all:
 input: expand(str(out/'{name}'),name=output_names)
rule compare:
 input: manifest=str(manifest),npz=paths,script=str(BASE/'render.py'),rscript=str(BASE/'statistics.R')
 output: expand(str(out/'{name}'),name=output_names)
 params: out=str(out),channel=config.get('channel','STAD'),sigma=float(config.get('sigma_bp',20)),permutations=int(config.get('permutations',99999)),seed=int(config.get('seed',20261110)),rscript=config.get('rscript','Rscript')
 shell: 'python3 {input.script:q} --manifest {input.manifest:q} --outdir {params.out:q} --channel {params.channel:q} --sigma-bp {params.sigma} --permutations {params.permutations} --seed {params.seed} --rscript {params.rscript:q}'
