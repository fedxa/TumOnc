import TumOnc as to
import argparse

parser = argparse.ArgumentParser(description='Run the analysis')
parser.add_argument("--mutfile", help="MAF with mutations",
                    default='BCC_MUT_combined_context.maf')
parser.add_argument("--covariates", help="gene covariates",
                    default='gene.covariates.txt')
parser.add_argument("--pval4", help="gene covariates",
                    default=False)
parser.add_argument("--minmuts", help="Minimal mutation count to analyse",
                    default=4)
parser.add_argument("--pmodel4", help="use pmodel 4",
                    default=0)
parser.add_argument("--pmodel_bagels", help="use bagels",
                    default=0)
parser.add_argument("--pmodel4bis", default=0)
parser.add_argument("--outfilebase", help="Basename for output file",
                    default="out")
parser.add_argument("--savepmodel", help="Just save the pmodel file",
                    default=False)
parser.add_argument("--loadpmodel",
                    help="Uses previously saved pmodel instead of calculation",
                    default=False)
parser.add_argument("--interactive", action="store_true")

args = parser.parse_args()
mutfile = args.mutfile
pval4 = args.pval4
if pval4:
    print("using critical statistics with 4 hits")
else:
    print("using critical statistics with 3 hits")

if args.pmodel4:
    print("Using pmodel4")
    pmod = to.pmodel.pmodel4(float(args.pmodel4))
elif args.pmodel_bagels:
    print("Using bagels as pmodel")
    pmod = to.pmodel.pmodel_bagels()
elif args.pmodel4bis:
    print("Using pmodel4bis avg=", args.pmodel4bis)
    pmod = to.pmodel.pmodel4bis(float(args.pmodel4bis))
else:
    print("Using pmodel3")
    pmod = to.pmodel.pmodel3()

if args.minmuts:
    minmuts = int(args.minmuts)
else:
    minmuts = 3
print("Using ", minmuts, " minmal mutations for analusis")

genes = to.read_genes_default(covariates=args.covariates)

print("Reading mutations from %s" % mutfile)
bcc_all=to.read_mutations(mutfile)

# Select only Missense mutations
mis_mut = to.clean_mutations(bcc_all[bcc_all['Variant_Classification']=='Missense_Mutation'],genes)
sil_mut = to.clean_mutations(bcc_all[bcc_all['Variant_Classification']=='Silent'],genes)
snp_mut = to.clean_mutations(bcc_all[bcc_all['Variant_Type']=='SNP'],genes)

allc5=to.context.get_ctx_bins(snp_mut)

if args.loadpmodel:
    print('Loading pmodel from %s instead of estimating' % args.loadpmodel)
    snp_mut_bin = to.binned_ns_load(args.loadpmodel)
else:
    snp_mut_bin = to.binned_ns(snp_mut, genes).bin(pmodel=pmod)

if args.savepmodel:
    print("Saving calculatted pmodel to "+args.savepmodel)
    snp_mut_bin.save(args.savepmodel)
    exit()

print('starting analysis')

mis_mut_bin = to.binned_ns(mis_mut, genes).bin(pmodel=0)
sil_mut_bin = to.binned_ns(sil_mut, genes).bin(pmodel=0)

mis_out_file = args.outfilebase+'_mis.csv'


if not args.interactive:
    mis_res, mis_hits = to.search_oncogenes(mis_mut_bin, snp_mut_bin,
                                            pval4=pval4, minmuts=minmuts)
    print("saving missense results to %s" % mis_out_file)
    mis_res.to_csv(mis_out_file, sep='\t')
    mis_res.head(20)

    sil_out_file = args.outfilebase+'_sil.csv'
    sil_res, sil_hits = to.search_oncogenes(sil_mut_bin, snp_mut_bin,
                                            pval4=pval4)
    print("saving silent results (control) to %s" % sil_out_file)
    sil_res.to_csv(sil_out_file, sep='\t')
    sil_res.head(20)
