from collections import defaultdict

try:
    from regex import regex
except:
    import regex
from config import *


def get_units(s, kmers):
    edges = defaultdict(int)
    vertices = defaultdict(int)

    g = defaultdict(list)
    for kmer, occ in kmers.items():
        v1 = kmer[:-1]
        v2 = kmer[1:]
        edges[(v1,v2)] += occ
        vertices[v1] += occ//2
        vertices[v2] += occ//2
        g[v1].append(v2)

    units = defaultdict(int)
    while max(edges.values()) > 0:
        kmers = [v for v,i in vertices.items() if i>0]
        if not kmers:
            break
        max_edge = max(edges.items(), key=lambda k: k[1])[0]
        start_v = min(kmers)
        unit = [start_v]
        v = start_v
        path_cov = max(edges.values())
        added_edges = []
        while g[v]:
            cur_str = unit[0] + ''.join([x[-1] for x in unit[1:]])
            next_f, next_v = 0, g[v][0]
            for i in g[v]:
                if vertices[i] > next_f and cur_str + i[-1] in s:
                    next_f,next_v=vertices[i],i
            path_cov = min(path_cov, edges[(v, next_v)])
            added_edges.append((v,next_v))
            v = next_v
            if next_v == start_v:
                break
            if next_f == 0:
                break
            unit.append(next_v)

        vertices[unit[0]] -= 1
        for i in range(1, len(unit)):
            edges[(unit[i - 1], unit[i])] -= 1
            if edges[(unit[i - 1], unit[i])] == 0:
                g[unit[i - 1]].remove(unit[i])
            vertices[unit[i]] -= 1
        edges[(unit[-1], unit[0])] -= 1
        if edges[(unit[-1], unit[0])] == 0:
            g[unit[-1]].remove(unit[0])
        if v == start_v:
            unit_str = unit[0] + ''.join([x[-1] for x in unit[1:-(TMER_SIZE-2)]])
            unit_start = unit_str.index(min(unit_str))
            unit_str = unit_str[unit_start:] + unit_str[:unit_start]
            s = s.replace(unit_str, '', 1)
            if len(unit) >= 3:
                units[unit_str] += 1
        else:
            unit_str = unit[0] + ''.join([x[-1] for x in unit[1:]])
            s = s.replace(unit_str, '', 1)
            if len(unit) >= 3:
                units[unit_str] += 1

    units_occ = [(unit,occ) for unit, occ in units.items()]
    units_occ.sort(key=lambda x:len(x[0]), reverse=True)
    return units_occ


def analyze_unit_structure(monomer_structure):
    gaps = []
    unit_structure = []
    monomer_structure.sort(key=lambda x:x[1])
    for i in range(1, len(monomer_structure)):
        if abs(monomer_structure[i][1] - monomer_structure[i-1][2]) > MONOMER_GAP_SIZE:
            gaps.append((i, min(monomer_structure[i-1][2], monomer_structure[i-1][1]), max(monomer_structure[i-1][2], monomer_structure[i-1][1]), monomer_structure[i-1][0]))
    monomers_str = ''.join([x[0] for x in monomer_structure])
    tmers = defaultdict(int)
    for i in range(len(monomers_str) - TMER_SIZE + 1):
        tmers[monomers_str[i:i + TMER_SIZE]] += 1
    if not tmers:
        return unit_structure

    units = get_units(monomers_str, tmers)
    if not units:
        return unit_structure
    for unit, occ in units:
        for unit_occ in regex.finditer(unit, monomers_str):
            unit_start, unit_end = unit_occ.span()
            unit_structure.append((monomer_structure[unit_start][1], monomer_structure[unit_end-1][2], unit))
            monomers_str = monomers_str[:unit_start]+'N'*(unit_end-unit_start)+monomers_str[unit_end:]
    unit_structure.sort(key=lambda x:x[0])
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
    return unit_structure


