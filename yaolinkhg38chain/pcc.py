from yaolinkhg38chain import hg38toyaochainlift
from yaolinkhg38chain import yaotohg38chainlift
from yaolinkhg38chain import chrtrans
from yaolinkhg38chain.arguments import get_args
import os


def main(arglist=None):
    # 1. 获取当前脚本文件的绝对路径（比如：/home/.../yaolinkhg38pcc/yaolinkhg38chain/hg38toyaochainlift.py）
    script_path = os.path.abspath(__file__)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    args = get_args(arglist)
    chrom = args.chr
    if args.pos:
        if args.parent=="mat":
            paf_file = f"{base_dir}/../data/yaomat2hg38v1.0.pcc"
            args.chr = chrtrans.yaomatchromosome_dict_rev[chrom]
        elif args.parent=="pat":
            paf_file = f"{base_dir}/../data/yaopat2hg38v1.0.pcc"
            args.chr = chrtrans.yaopatchromosome_dict_rev[chrom]
        if args.src == 'hg38' and args.dest == 'yao':
            args.chr = chrtrans.hg38chromosome_dict_rev[chrom]
            hg38toyaochainlift.hg38toyaochainlift(paf_file, args.chr, args.pos,args.parent)
        elif args.src == 'yao' and args.dest == 'hg38':
            yaotohg38chainlift.yaotohg38chainlift(paf_file, args.chr, args.pos,args.parent)
    elif args.Seg:
        if args.parent=="mat":
            paf_file = f"{base_dir}/../data/yaomat2hg38v1.0.pcc"
            args.chr = chrtrans.yaomatchromosome_dict_rev[chrom]
        elif args.parent=="pat":
            paf_file = f"{base_dir}/../data/yaopat2hg38v1.0.pcc"
            args.chr = chrtrans.yaopatchromosome_dict_rev[chrom]
        if args.src == 'hg38' and args.dest == 'yao':
            args.chr = chrtrans.hg38chromosome_dict_rev[chrom]
            hg38toyaochainlift.hg38toyaochainlift_locus(paf_file, args.chr, args.Seg,args.parent)
        elif args.src == 'yao' and args.dest == 'hg38':
            yaotohg38chainlift.yaotohg38chainlift_locus(paf_file, args.chr, args.Seg,args.parent)
    elif args.posfile:
        if args.parent=="mat":
            paf_file = f"{base_dir}/../data/yaomat2hg38v1.0.pcc"
            args.chr = chrtrans.yaomatchromosome_dict_rev[chrom]
        elif args.parent=="pat":
            paf_file = f"{base_dir}/../data/yaopat2hg38v1.0.pcc"
            args.chr = chrtrans.yaopatchromosome_dict_rev[chrom]
        if args.src == 'hg38' and args.dest == 'yao':
            args.chr = chrtrans.hg38chromosome_dict_rev[chrom]
            with open(args.posfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        chr, pos = line.split()
                        pos = int(pos)
                        hg38toyaochainlift.hg38toyaochainlift(paf_file, args.chr, pos,args.parent)
        elif args.src == 'yao' and args.dest == 'hg38':
            with open(args.posfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        chr, pos = line.split()
                        pos = int(pos)
                        yaotohg38chainlift.yaotohg38chainlift(paf_file, args.chr, pos,args.parent)
    elif args.Segfile:
        if args.parent=="mat":
            paf_file = f"{base_dir}/../data/yaomat2hg38v1.0.pcc"
            args.chr = chrtrans.yaomatchromosome_dict_rev[chrom]
        elif args.parent=="pat":
            paf_file = f"{base_dir}/../data/yaopat2hg38v1.0.pcc"
            args.chr = chrtrans.yaopatchromosome_dict_rev[chrom]
        if args.src == 'hg38' and args.dest == 'yao':
            args.chr = chrtrans.hg38chromosome_dict_rev[chrom]
            with open(args.Segfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        chr, start, end = line.split()
                        start = int(start)
                        end = int(end)
                        hg38toyaochainlift.hg38toyaochainlift_locus(paf_file, args.chr, (start, end),args.parent)
        elif args.src == 'yao' and args.dest == 'hg38':
            with open(args.Segfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        chr, start, end = line.split()
                        start = int(start)
                        end = int(end)
                        yaotohg38chainlift.yaotohg38chainlift_locus(paf_file, args.chr, (start, end),args.parent)
    

if __name__ == '__main__':
    main()

#def chrtrans(,chr):
    
