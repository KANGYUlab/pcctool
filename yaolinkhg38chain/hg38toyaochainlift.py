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
    """
    ops = []
    parts = re.findall(r'(\d+)([MID])', cg_str)
    for length_str, op in parts:
        ops.append((op, int(length_str)))
    return ops
def ref_to_query_pos_via_cg(ref_pos, rstart, rend, qstart, qend, strand, cg_str):
    # if not (rstart <= ref_pos < rend):
    #     return None

    ops = parse_cigar(cg_str)
    r_curr = rstart
    q_curr = qstart  # q_curr 始终代表正向比对的进度

    # --- 提前根据 strand 确定转换逻辑 ---
    # 我们计算一个相对于 qstart 的位移 dist = q_pos_forward - qstart
    # 如果是正链: final_q = qstart + dist
    # 如果是负链: final_q = (qend - 1) - dist
    is_forward = (strand == '+')

    for op, length in ops:
        if op == 'M':
            # 边界判定：使用左闭右开 [r_curr, r_curr + length)
            if r_curr <= ref_pos < r_curr + length:
                offset = ref_pos - r_curr
                dist = (q_curr - qstart) + offset # 记录在 query 比对区间内走了多远
                
                if is_forward:
                    return qstart + dist
                else:
                    # 负链：从 qend-1 开始倒扣
                    # 对应你之前图片的公式：qend - (q_pos_forward - qstart) - 1
                    return (qend - 1) - dist
            
            r_curr += length
            q_curr += length

        elif op == 'D':
            if r_curr <= ref_pos < r_curr + length:
                return None  # 落在缺失区
            r_curr += length

        elif op == 'I':
            # Insertion 只增加 query 进度，不增加 ref 进度
            q_curr += length
            
        elif op in ('S', 'H', 'N'):
            # 处理可能存在的其他标签
            if op != 'N': q_curr += length
            if op != 'I': r_curr += length
            
        # 性能优化：如果已经超过了目标位置，提前结束
        if r_curr > ref_pos:
            break

    return None

def parse_paf_with_cg(line):
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None
    qname = fields[0]
    qstart = int(fields[2])
    qend = int(fields[3])
    strand = fields[4]
    tname = fields[5]
    rstart = int(fields[7])
    rend = int(fields[8])
    confidence = float(fields[15])
    # 搜索 cg:Z: 字段
    cg_tag = None
    for field in fields[12:]:
        if field.startswith('cg:Z:'):
            cg_tag = field[5:]
            break
    if cg_tag is None:
        return None

    return {
        'qname': qname,
        'qstart': qstart,
        'qend': qend,
        'strand': strand,
        'tname': tname,
        'rstart': rstart,
        'rend': rend,
        'cg': cg_tag,
        'confidence': confidence
    }

def hg38toyaochainlift(paf_file, ref_chrom, ref_pos_0based,parent):
    found = False
    with open(paf_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            aln = parse_paf_with_cg(line)
            if not aln:
                continue
            if aln['tname'] != ref_chrom:
                continue
            if not (aln['rstart'] <= ref_pos_0based < aln['rend']):
                continue

            qpos = ref_to_query_pos_via_cg(
                ref_pos=ref_pos_0based,
                rstart=aln['rstart'],
                rend=aln['rend'],
                qstart=aln['qstart'],
                qend=aln['qend'],
                strand=aln['strand'],
                cg_str=aln['cg']
            )
            if aln['strand'] == '+':
                qpos = qpos
            else:
                qpos = qpos + 1
            if qpos is not None:
                print(f"hg38 {chrtrans.hg38chromosome_dict[ref_chrom]}:{ref_pos_0based} → yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[aln['qname']]}:{qpos} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
                found = True
                # 可选：只返回第一个匹配，或继续找多个
                break  # 如果允许多次比对，可注释掉 break

    if not found:
        print(f"No alignment with 'cg' tag covers {chrtrans.hg38chromosome_dict[ref_chrom]}:{ref_pos_0based}")
# if __name__ == "__main__":
#     # if len(sys.argv) != 4:
#     #     print("Usage: python paf_ref2query_cg.py <paf_file> <ref_chrom> <ref_position_0based>")
#     #     print("Note: PAF must be generated with minimap2 --cg")
#     #     sys.exit(1)

#     paf_path = "../data/yaomat20k2hg38blockendv1.txt"
#     ref_chrom = "NC_000002.12"
#     ref_pos = int(110160255)
def hg38toyaochainlift_locus(paf_file, chrom, locus,parent):
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
                if aln['tname'] != chrom:
                    continue
                if not (aln['rstart'] <= ref_pos_0based < aln['rend']):
                    continue

                qpos = ref_to_query_pos_via_cg(
                    ref_pos=ref_pos_0based,
                    rstart=aln['rstart'],
                    rend=aln['rend'],
                    qstart=aln['qstart'],
                    qend=aln['qend'],
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
            print(f"No alignment with 'cg' tag covers {getattr(chrtrans, f"yao{parent}chromosome_dict")[aln['qname']]}:{ref_pos_0based}")
            liftlocus.append(None)
    if aln['strand'] == '+':
        print(f"hg38 {chrtrans.hg38chromosome_dict[chrom]}:{locus} → yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[aln['qname']]}:{liftlocus} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
        return liftlocus
    else:
        print(f"hg38 {chrtrans.hg38chromosome_dict[chrom]}:{locus} → yao{parent} {getattr(chrtrans, f"yao{parent}chromosome_dict")[aln['qname']]}:{[liftlocus[1]+1,liftlocus[0]+1]} (strand={aln['strand']}, ChainQ={aln['confidence']:.4f})")
        return [liftlocus[1]+1,liftlocus[0]+1]
