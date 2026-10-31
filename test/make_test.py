import csv, random, gzip, sys
random.seed(1)
N=40; S=[f"S{i:03d}" for i in range(N)]
rows=list(csv.DictReader(open("resources/app_variants_v2.tsv"),delimiter="\t"))
recs={}
for r in rows:
    if r["pos38"]: recs[(r["chrom"],int(r["pos38"]),r["ref"],r["alt"])]=r["rsid"]
lines=[]
flipped="rs1057141"   # TAP1 I333V written with REF/ALT swapped
for (c,p,ref,alt),rs in sorted(recs.items(),key=lambda x:(int(x[0][0]),x[0][1])):
    if rs==flipped: ref,alt=alt,ref
    af=random.uniform(0.05,0.5)
    r2 = 0.1 if rs=="rs17587" else 0.95          # PSMB9 R60H fails R2
    gts=[]
    for i,s in enumerate(S):
        a=int(random.random()<af); b=int(random.random()<af)
        g=f"{a}|{b}"
        if rs=="rs2228396" and i==3: g=".|."    # missing call in TAP2
        gts.append(g)
    lines.append([f"chr{c}",str(p),rs,ref,alt,".","PASS",f"R2={r2}","GT"]+gts)
# extra non-catalogue variant that must NOT be extracted
lines.append(["chr6","32850000","rsNOTME","A","G",".","PASS","R2=0.99","GT"]+["0|1"]*N)
lines.sort(key=lambda l:(int(l[0][3:]),int(l[1])))
hdr=["##fileformat=VCFv4.2",'##INFO=<ID=R2,Number=1,Type=Float,Description="r2">','##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">']
hdr+=[f"##contig=<ID=chr{c}>" for c in [1,5,6,9,12,14,15,16,17,19]]
hdr.append("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+S))
open("test/wgs.vcf","w").write("\n".join(hdr+["\t".join(l) for l in lines])+"\n")
# HLA VCF
h=["##fileformat=VCFv4.2",'##INFO=<ID=DR2,Number=1,Type=Float,Description="dr2">','##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',"##contig=<ID=6>",
   "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+S)]
for k,(nm,pos) in enumerate([("HLA_A*01",29942000),("HLA_A*01:01",29942001),("HLA_B*07",31353000),("HLA_B*07:02",31353001)]):
    h.append("\t".join(["6",str(pos),nm,"A","T",".","PASS","DR2=0.95","GT"]+[random.choice(["0|0","0|1","1|1","1|0"]) for _ in S]))
open("test/hla.vcf","w").write("\n".join(h)+"\n")
# fake hlafactor
fake=r'''#!/usr/bin/env python3
import sys,os,subprocess
a=sys.argv; i=a[a.index("-i")+1]; o=a[a.index("-o")+1]; m=a[-1][1:]
samples=subprocess.run(["bcftools","query","-l",i],capture_output=True,text=True).stdout.split()
d=os.path.dirname(o)
def w(fn,cols):
    with open(os.path.join(d,m+fn),"w") as f:
        f.write("\t".join(["SampleID"]+cols)+"\n")
        for s in samples: f.write("\t".join([s]+["1"]*len(cols))+"\n")
{"het":lambda:w("_heterozygosity.txt",["A","B","C"]),
 "hetf":lambda:w("_functional_zygosity.txt",["HLA_A_FH","HLA_B_FH","HLA_C_FH"]),
 "exp":lambda:(w("_exp_A.txt",["expression_value"]),w("_exp_C.txt",["expression_value"])),
 "sup":lambda:w("_supertypes.txt",["A01","B07"]),
 "tapasin":lambda:[w(f"_taps_{x}.txt",[f"taps{x}_value" if x!="global" else "taps_global_value"]) for x in ["A","B","C","global"]],
 "lignk":lambda:w("_HLA_variant_NK.txt",["HLA_C_ligand","HLA_B_ligand"])}[m]()
'''
open("test/fake_hlafactor","w").write(fake)
