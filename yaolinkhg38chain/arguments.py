import argparse
import sys

def get_args(arglist):
    parser = argparse.ArgumentParser(description='pcctool: hg38 <-> yao mapping')
    
    # 1. 转换方向 (从哪到哪)
    parser.add_argument('-src', choices=['hg38', 'yao'], required=True, help='源基因组')
    parser.add_argument('-dest', choices=['hg38', 'yao'], required=True, help='目标基因组')
    
    # 2. 亲本判断 (只有 yao 需要判断亲本)
    parser.add_argument('-p', '--parent', choices=['mat', 'pat'], required=True, help='亲本: mat(母本) 或 pat(父本)')
    
    # 3. 坐标信息
    parser.add_argument('-chr', required=True, help='染色体[chr1,chr2,...]')
    parser.add_argument('-pos', type=int, required=False, help='位置')
    parser.add_argument('-Seg', nargs=2, type=int, required=False, metavar=('START', 'END'), 
                        help='Segment区域范围（两个整数表示区间，如 100 200）')
    # 批量位置文件：每行格式预计为 "chr pos"
    parser.add_argument('-posfile', type=str, required=False, help='包含坐标列表的TXT文件路径')

    # 批量区间文件：每行格式预计为 "chr start end"
    parser.add_argument('-Segfile', type=str, required=False, help='包含区间列表的TXT文件路径')
    parser.add_argument('-V', '--version', help='show program version', action='version', version='1.0')
    args = parser.parse_args(arglist)
    if args.src == args.dest:
        sys.exit("错误: 源和目标不能相同")
    return parser.parse_args()
