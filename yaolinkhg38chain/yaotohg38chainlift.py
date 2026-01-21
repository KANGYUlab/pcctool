#!/usr/bin/env python3
import sys
import re
from yaolinkhg38chain import chrtrans
# ==========================================================
# 辅助函数 - CIGAR 
# ==========================================================

def parse_cigar(cg_str):
    """
    解析 CIGAR 字符串（如 "4M2D10M2I"）为操作列表。
    返回: [('M', 4), ('D', 2), ('M', 10), ('I', 2), ...]
    注意: 仅解析 M, I, D 操作，忽略其他（如 N, S, H 等）。
    """
    ops = []
    # 扩展匹配，以支持 Minimap2 cg:Z: 中可能出现的 N (=, X) 等，但我们只处理 MID
    parts = re.findall(r'(\d+)([MIDSHN=X])', cg_str)
    for length_str, op in parts:
        # 仅将 MID 计入坐标推进
        ops.append((op, int(length_str)))
    return ops

def parse_paf_with_cg(line):
    """从 PAF 行解析出关键字段和 cg:Z: 标签"""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None
        
    try:
        data = {
            'qname': fields[0],
            'qstart': int(fields[2]),
            'qend': int(fields[3]),
            'strand': fields[4],
            'tname': fields[5],
            'rstart': int(fields[7]),
            'rend': int(fields[8]),
            'confidence': float(fields[15]),
        }
    except ValueError:
        # 字段非数字，跳过
        return None

    # 搜索 cg:Z: 字段
    data['cg'] = None
    for field in fields[12:]:
        if field.startswith('cg:Z:'):
            data['cg'] = field[5:]
            break
            
    if data['cg'] is None:
        return None
    return data

# ==========================================================
# 核心函数：Query 位置推导 Reference 位置
# ==========================================================
#Query->ref
def query_to_ref_pos_via_cg(query_pos, qstart, qend, rstart, rend, strand, cg_str):


    # 1. 处理负链坐标转换
    if strand == '+':
        target_q = query_pos
    else:
        # 负链逻辑：将递减的 Query 坐标映射为随 CIGAR 步进递增的等效坐标
        # 公式：qstart + (qend - 1 - query_pos)
        target_q = qstart + (qend - 1 - query_pos)
    # 初始化进度
    ops = parse_cigar(cg_str)
    q_curr = qstart
    r_curr = rstart
    for op, length in ops:
        if op == 'M' or op == '==' or op == 'X':
            # Match: ref 和 query 同时前进
            # 建议：这里使用 q_curr <= target_q < q_curr + length 更稳健
            if q_curr <= target_q <= q_curr + length:
                offset = target_q - q_curr
                ref_pos_forward = r_curr + offset
                return ref_pos_forward
            
            q_curr += length
            r_curr += length
        elif op == 'I':
            # Insertion in Query: query 有, reference 无
            if q_curr <= target_q <= q_curr + length:
                # Query 位置落在 insertion 区 -> 无对应 Reference 位置
                return None
            q_curr += length
            # r_curr 不变
        elif op == 'D':
            # Deletion in Query: reference 有, query 无
            r_curr += length
            # q_curr 不变
        elif op in ('S', 'H'):
            # 软/硬裁剪：通常只影响 query 坐标
            q_curr += length
            
        elif op == 'N':
            # 跳跃（内含子）
            r_curr += length
    return None
# ==========================================================
# 主逻辑 (仿照您的 main 函数结构)
# ==========================================================

def yaotohg38chainlift(paf_file, query_chrom, query_pos_0based, parent):
    found = False
    
    # 

    with open(paf_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
                
            aln = parse_paf_with_cg(line)
            if not aln:
                continue
                
            # 1. 筛选 Query 染色体
            if aln['qname'] != query_chrom:
                continue
                
            # 2. 筛选位置所在的比对块
            if not (aln['qstart'] <= query_pos_0based < aln['qend']):
                continue

            # 3. 进行坐标转换
            ref_pos = query_to_ref_pos_via_cg(
                query_pos=query_pos_0based,
                qstart=aln['qstart'],
                qend=aln['qend'],
                rstart=aln['rstart'],
                rend=aln['rend'],
                strand=aln['strand'],
                cg_str=aln['cg']
            )
            if aln['strand'] == "+":
                ref_pos=ref_pos
            else:
                ref_pos=ref_pos+1
            if ref_pos is not None:
                print(f"yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[query_chrom]}:{query_pos_0based} → hg38 {chrtrans.hg38chromosome_dict[aln['tname']]}:{ref_pos} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
                found = True
                # 可选：只返回第一个匹配，或继续找多个
                break 

    if not found:
        print(f"No alignment with 'cg' tag covers yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[query_chrom]}:{query_pos_0based}")

# if __name__ == "__main__":
#     # 请根据您的实际文件和位置替换以下示例值
#     #paf_path = "../data/yaomat20k2hg38blockendv1.txt"
#     paf_path = "../data/yaomat20k2hg38blockendv1.txt"
#     query_chrom = "chr2_hap1|110522120|110611322|block375" # 假设 Query 序列名为 chr1_hap1Query chr1_hap1|71461|227856|block1:71462 (strand=+)
#     query_pos = 110574412     # 假设要转换的 0-based 位置
#     #chr1_hap1|71461|227856|block1:121517
#     #chr2_hap1|110522120|110611322|block375:110574412 (strand=-)
#     # 命令行参数用法（如果您需要）
#     # if len(sys.argv) == 4:
#     #     paf_path = sys.argv[1]
#     #     query_chrom = sys.argv[2]
#     #     try:
#     #         query_pos = int(sys.argv[3])
#     #     except ValueError:
#     #         print("Error: query_position must be an integer.")
#     #         sys.exit(1)

#     main_query_to_ref_example(paf_path, query_chrom, query_pos)
def yaotohg38chainlift_locus(paf_file, chrom, locus, parent):
    liftlocus=[]
    for position in [locus[0], locus[1]]:
        ref_pos_0based = position
        found = False
        with open(paf_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                aln = parse_paf_with_cg(line)
                if not aln:
                    continue
                if aln['qname'] != chrom:
                    continue
                if not (aln['qstart'] <= ref_pos_0based < aln['qend']):
                    continue

                qpos = query_to_ref_pos_via_cg(
                    query_pos=position,
                    qstart=aln['qstart'],
                    qend=aln['qend'],
                    rstart=aln['rstart'],
                    rend=aln['rend'],
                    strand=aln['strand'],
                    cg_str=aln['cg']
                )
                if qpos is not None:
                    #print(f"Reference {chrom}:{ref_pos_0based} → Query {aln['qname']}:{qpos} (strand={aln['strand']})")                    
                    found = True
                    liftlocus.append(qpos)
                    # 可选：只返回第一个匹配，或继续找多个
                    break  # 如果允许多次比对，可注释掉 break

        if not found:
            print(f"No alignment with 'cg' tag covers yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[chrom]}:{ref_pos_0based}")
            liftlocus.append(None)
    if aln['strand'] == '+':
        print(f"yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[chrom]}:{locus} → hg38 {chrtrans.hg38chromosome_dict[aln['tname']]}:{liftlocus} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
        return liftlocus
    else:
        print(f"yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[chrom]}:{locus} → hg38 {chrtrans.hg38chromosome_dict[aln['tname']]}:{[liftlocus[1]+1,liftlocus[0]+1]} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
        return [liftlocus[1]+1,liftlocus[0]+1]