import re
from collections import defaultdict

import numpy as np
try:
    from regex import regex
except:
    import regex
from config import *


def analyze_unit_structure(monomer_structure):
    gaps = []
    unit_structure = []
    monomer_structure.sort(key=lambda x:x[1])
    for i in range(1, len(monomer_structure)):
        if abs(monomer_structure[i][1] - monomer_structure[i-1][2]) > MONOMER_GAP_SIZE:
            gaps.append((i, min(monomer_structure[i-1][2], monomer_structure[i-1][1]), max(monomer_structure[i-1][2], monomer_structure[i-1][1]), monomer_structure[i-1][0]))
    monomers_str = ''.join([x[0] for x in monomer_structure])
    tmers = defaultdict(list)
    for i in range(len(monomers_str) - TMER_SIZE + 1):
        tmers[monomers_str[i:i + TMER_SIZE]].append(i)
    mean_occ = np.percentile([len(occ) for t, occ in tmers.items()], 95)
    distances = []
    for tmer, occ in tmers.items():
        if len(occ) < mean_occ:
            continue
        dist = [occ[i] - occ[i - 1] for i in range(1, len(occ))]
        if dist:
            distances.append(np.median(dist))

    if not distances:
        return unit_structure, None

    d = int(np.median(distances))
    all_monomers = defaultdict(int)
    for i in range(len(monomers_str) - d + 1):
        all_monomers[monomers_str[i:i+d]] += 1

    monomers_pattern = max(all_monomers, key=lambda k: all_monomers[k])
    if 'A' in monomers_pattern:
        monomers_pattern = monomers_pattern[monomers_pattern.index('A'):] + monomers_pattern[:monomers_pattern.index('A')]
    #print(monomers_pattern)
    #monomers_pattern = "ABCDEFGHIJKL"
    prev_start, prev_end = 0, 0
    unit_end = 0
    universal_pattern = r'(' + re.escape(monomers_pattern) + r'|' + re.escape(monomers_pattern[::-1]) + r'){e<=0}'
    for m in regex.finditer(universal_pattern, monomers_str, regex.ENHANCEMATCH):
        def reverse_substr(substr, start, end):
            if start == 0:
                return substr[end::-1]
            else:
                return substr[end:start-1:-1]
        unit_start, unit_end = m.span()
        if unit_start - prev_end > 2:
            if unit_start - prev_end <= 5:
                unit_structure.append((monomer_structure[prev_end][1], monomer_structure[unit_start - 1][2],
                                       monomers_str[prev_end:unit_start]))
            else:
                wrong_substr = monomers_str[prev_end:unit_start]
                str_len, i = 1, 0
                while i < len(wrong_substr) - 1:
                    while i + str_len < len(wrong_substr) and \
                            (wrong_substr[i:i + str_len] in monomers_pattern or reverse_substr(wrong_substr, i, i + str_len) in monomers_pattern):
                        str_len += 1
                    if wrong_substr[i:i + str_len] not in monomers_pattern and reverse_substr(wrong_substr, i, i + str_len) not in monomers_pattern:
                        str_len -= 1
                    if str_len >= 2:
                        unit_structure.append((monomer_structure[prev_end + i][1],
                                               monomer_structure[prev_end + i + str_len - 1][2],
                                               monomers_str[prev_end + i:prev_end + i + str_len]))
                    i = i + max(1, str_len)
                    str_len = 1
                if str_len >= 2:
                    unit_structure.append((monomer_structure[prev_end + i][1],
                                           monomer_structure[prev_end + i + str_len - 1][2],
                                            monomers_str[prev_end + i:prev_end + i + str_len]))
        elif unit_structure:
            unit = unit_structure[-1]
            unit_structure[-1] = (unit[0], monomer_structure[unit_start - 1][2], unit[2] + monomers_str[prev_end:unit_start])
        unit_structure.append((monomer_structure[unit_start][1], monomer_structure[unit_end - 1][2],
                                   monomers_str[unit_start:unit_end]))
        prev_start, prev_end = unit_start, unit_end
    if unit_end + 2 < len(monomer_structure):
        unit_structure.append((monomer_structure[unit_end + 1][1], monomer_structure[-1][2],
                                monomers_str[unit_end + 1:]))

    prev_s, prev_e = 0, 0
    i = 0
    gaps = sorted(list(gaps))
    gaps_n = 0
    while i < len(unit_structure):
        unit_start, unit_end, hor = unit_structure[i]
        while gaps_n < len(gaps) and unit_start > gaps[gaps_n][1]:
            gaps_n += 1
        if gaps_n < len(gaps) and unit_start < gaps[gaps_n][1] and gaps[gaps_n][2] < unit_end:
            #print(prev_e,unit_start,i,len(unit_structure),hor, gaps)
            #gaps.append((prev_e, unit_start))
            unit_structure.remove(unit_structure[i])
            if len(hor[:hor.index(gaps[gaps_n][3]) + 1]) > 1:
                unit_structure.insert(i, (unit_start, gaps[gaps_n][1], hor[:hor.index(gaps[gaps_n][3]) + 1]))
                i += 1
            if len(hor[hor.index(gaps[gaps_n][3]) + 1:]) > 1:
                unit_structure.insert(i + 1, (gaps[gaps_n][2], unit_end, hor[hor.index(gaps[gaps_n][3]) + 1:]))
                i += 1
        elif unit_start - prev_e > MONOMER_GAP_SIZE:
            pass
            #print("GAP", prev_e, unit_start)
            #gaps.append((prev_e, unit_start))
        else:
            i += 1
        prev_s, prev_e = unit_start, unit_end
    return unit_structure, monomers_pattern


