#!/usr/bin/env python3
"""Run from any directory. Compiler flags, seeds, and raw measurements are retained."""
import subprocess,pathlib,csv,tempfile,sys,statistics
base=pathlib.Path(__file__).resolve().parent
out=pathlib.Path(sys.argv[1] if len(sys.argv)>1 else base/'benchmark.csv')
n=int(sys.argv[2]) if len(sys.argv)>2 else 200000
repeats=int(sys.argv[3]) if len(sys.argv)>3 else 3
with tempfile.TemporaryDirectory() as tmp, out.open('w') as f:
    writer=csv.writer(f);writer.writerow(['implementation','monoid','repeat','n','mode','build_s','operations_s','peak_rss_kib','checksum','build_cpu_s','operations_cpu_s'])
    for affine in [False,True]:
        binaries={}
        for old in [True,False]:
            exe=pathlib.Path(tmp)/f'bench_{affine}_{old}'
            cmd=['g++','-std=c++20','-O2','-DNDEBUG']+(['-DOLD'] if old else [])+(['-DAFFINE'] if affine else [])+[str(base/'bench.cpp'),'-o',str(exe)]
            subprocess.run(cmd,check=True);binaries[old]=exe
        for mode in range(6):
            for repeat in range(repeats):
                for old in ([True,False] if repeat%2==0 else [False,True]):
                    if old and mode in [2,3]:continue
                    result=subprocess.check_output([str(binaries[old]),str(n),str(mode)],text=True).strip().split(',')
                    writer.writerow(['old' if old else 'new','affine' if affine else 'sum',repeat]+result);f.flush()
print(out)
