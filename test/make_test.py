import csv, random, gzip, sys
random.seed(1)
N=int(__import__("os").environ.get("NS","40")); S=[f"S{i+1}" for i in range(N)] if __import__("os").environ.get("REALNAMES") else [f"S{i:03d}" for i in range(N)]
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
# Stand-in for hlafactor (https://github.com/QiJingS/hlafactor) that writes files
# with the same names and formats as hlafactor v0.0.0, including text genotypes
# ("Het"/"Hom", "C1/C2") and "nan" missing values.
import sys,os,subprocess,random
a=sys.argv; i=a[a.index("-i")+1]; o=a[a.index("-o")+1]; m=a[-1][1:]
S=subprocess.run(["bcftools","query","-l",i],capture_output=True,text=True).stdout.split()
d=os.path.dirname(o) or "."; random.seed(m)
def w(fn,cols,f):
    with open(os.path.join(d,m+fn),"w") as fh:
        fh.write("\t".join(cols)+"\n")
        for k,s in enumerate(S): fh.write("\t".join([s]+f(k))+"\n")
miss=lambda k,x: "nan" if k==0 else x
if m=="het": w("_heterozygosity.txt",["SampleID","HLA_A","HLA_B","HLA_C",""],lambda k:[random.choice(["Het","Hom"]),"NA" if k==0 else "Het","Het"])
elif m=="hetf": w("_functional_zygosity.txt",["sampleid","HLA_A_allele","HLA_A_FH","HLA_B_allele","HLA_B_FH","HLA_C_allele","HLA_C_FH"],lambda k:["x","0.5","x",miss(k,"0.6"),"x","0.7"])
elif m=="exp":
    w("_exp_A.txt",["SampleID","locus","alleles","miss_match","expression_value"],lambda k:["HLA_A","x","NA","-0.2"])
    w("_exp_C.txt",["SampleID","locus","alleles","miss_match","expression_value"],lambda k:["HLA_C","x","NA","300"])
elif m=="sup": w("_supertypes.txt",["SampleID","A01","A02","B07",""],lambda k:["1","1","0"])
elif m=="tapasin":
    for x in "ABC": w(f"_taps_{x}.txt",["SampleID","locus","alleles","miss_match",f"taps{x}_value"],lambda k:["HLA_"+x,"x","NA",miss(k,"0.8")])
    w("_taps_global.txt",["SampleID","tapsA_value","tapsB_value","tapsC_value","taps_global_value"],lambda k:["0.8","0.8","0.8",miss(k,"1.2")])
elif m=="lignk": w("_HLA_variant_NK.txt",["SampleID","HLA_C_ligand","HLA_B_ligand"],lambda k:[random.choice(["C1/C1","C1/C2","C2/C2"]),"BW6/." if k==0 else random.choice(["BW4/BW6","BW6/BW6"])])
'''
open("test/fake_hlafactor","w").write(fake)
