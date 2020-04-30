#!/usr/bin/env python
import argparse


def delta2paf():
    """
    Essentially a translation of paftools.js delta2paf

    https://github.com/lh3/minimap2/tree/master/misc
    """
    parser = argparse.ArgumentParser(description="Convert a Nucmer delta file to a PAF file.")
    parser.add_argument("delta_file", metavar="<alns.delta>", type=str, help="delta file to convert.")

    args = parser.parse_args()
    delta_file = args.delta_file

    seen_gt = False
    with open(delta_file, "r") as f:
        for line in f:
            t = line.rstrip().split(" ")

            # Check to see if we have started a new ref/query pair
            if line.startswith(">"):
                rname, qname, rlen, qlen = t[0][1:], t[1], int(t[2]), int(t[3])
                seen_gt = True
                continue

            # Continue if we haven't reached alignments yet
            if not seen_gt:
                continue

            if len(t) == 7:
                for i in range(5):
                    t[i] = int(t[i])

                strand = True if (t[0] < t[1] and t[2] < t[3]) or (t[0] > t[1] and t[2] > t[3]) else False

                # Reference/query start/end
                rs = (t[0] - 1) if t[0] < t[1] else (t[1] - 1)
                re = t[1] if t[1] > t[0] else t[0]
                qs = (t[2] - 1) if t[2] < t[3] else (t[3] - 1)
                qe = t[3] if t[3] > t[2] else t[2]
                x, y = 0, 0
                NM = t[4]
                cigar = []

            elif len(t) == 1:
                d = int(t[0])
                if d == 0:
                    blen, cigar_str = 0, []
                    if re - rs - x != qe - qs - y:
                        print(re, rs, x)
                        print(qe, qs, y)
                        raise RuntimeError("inconsistent alignment")

                    cigar.append((re - rs - x) << 4)
                    for i in cigar:
                        blen += i >> 4
                        cigar_str.append(str((i >> 4)) + "MID"[i & 0xf])
                    out_strand = "+" if strand > 0 else "-"
                    print("\t".join([
                        qname,
                        str(qlen),
                        str(qs),
                        str(qe),
                        out_strand,
                        rname,
                        str(rlen),
                        str(rs),
                        str(re),
                        str(blen - NM),
                        str(blen),
                        "0",
                        "NM:i:" + str(NM),
                        "cg:Z:" + "".join(cigar_str)
                    ]))
                elif d > 0:
                    l = d - 1
                    x += l + 1
                    y += l
                    if l:
                        cigar.append(l << 4)
                    if len(cigar) > 0 and (cigar[-1] & 0xf) == 2:
                        cigar[-1] += 1 << 4
                    else:
                        cigar.append(1 << 4 | 2)  # deletion
                else:
                    l = -d - 1
                    x += l
                    y += l + 1
                    if l:
                        cigar.append(l << 4)
                    if len(cigar) > 0 and (cigar[-1] & 0xf) == 1:
                        cigar[-1] += 1 << 4
                    else:
                        cigar.append(1 << 4 | 1)  # insertion


if __name__ == "__main__":
    delta2paf()
