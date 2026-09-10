#!/usr/bin/env python3
"""Structural checks on a --homozyg group pool report.

The report repeats enough information to be checked against itself: every
pool's consensus region has to be contained in each member's run and the union
has to contain it, the CON row's count has to be the member count, and the
first group's size has to be its reference's NSIM plus one, since that group
is formed before anything has been assigned.
"""
import sys


def main():
    overlap_path, hom_path, want_pool_ct = sys.argv[1], sys.argv[2], int(sys.argv[3])
    runs = {}
    with open(hom_path) as f:
        header = f.readline().split()
        idx = {name.lstrip('#'): i for i, name in enumerate(header)}
        for line in f:
            g = line.split()
            runs.setdefault(g[idx['IID']], []).append(
                (int(g[idx['POS1']]), int(g[idx['POS2']])))
    pools = {}
    order = []
    with open(overlap_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        col = {name.lstrip('#'): i for i, name in enumerate(header)}
        for name in ('POOL', 'IID', 'POS1', 'POS2', 'NSIM', 'GRP'):
            if name not in col:
                sys.exit('missing %s column in the pool report' % name)
        # CON/UNION rows put the member count where IID normally goes, and the
        # marker itself in the first column after POOL.
        tag_col = col.get('FID', col['IID'])
        for line in f:
            g = line.rstrip('\n').split('\t')
            pool = g[col['POOL']]
            if pool not in pools:
                pools[pool] = {'members': [], 'con': None, 'union': None}
                order.append(pool)
            entry = (int(g[col['IID']]) if g[tag_col] in ('CON', 'UNION') else 0,
                     int(g[col['POS1']]), int(g[col['POS2']]))
            if g[tag_col] == 'CON':
                pools[pool]['con'] = entry
            elif g[tag_col] == 'UNION':
                pools[pool]['union'] = entry
            else:
                pools[pool]['members'].append(
                    (g[col['IID']], int(g[col['POS1']]), int(g[col['POS2']]),
                     int(g[col['NSIM']]), g[col['GRP']]))
    if len(pools) != want_pool_ct:
        sys.exit('%d pools, expected %d' % (len(pools), want_pool_ct))
    if order != sorted(order, key=lambda p: int(p[1:])):
        sys.exit('pool IDs are not in increasing order')
    prev_size = None
    for pool in order:
        info = pools[pool]
        members = info['members']
        if info['con'] is None or info['union'] is None:
            sys.exit('%s: missing CON or UNION row' % pool)
        con_ct, con_start, con_end = info['con']
        union_ct, union_start, union_end = info['union']
        if con_ct != len(members) or union_ct != len(members):
            sys.exit('%s: CON/UNION member count %d/%d, saw %d rows' %
                     (pool, con_ct, union_ct, len(members)))
        if prev_size is not None and len(members) > prev_size:
            sys.exit('%s: pools are not in decreasing size order' % pool)
        prev_size = len(members)
        for iid, pos1, pos2, nsim, grp in members:
            if pos1 > con_start or pos2 < con_end:
                sys.exit('%s: %s does not cover the consensus region' % (pool, iid))
            if pos1 < union_start or pos2 > union_end:
                sys.exit('%s: %s lies outside the union region' % (pool, iid))
            if (pos1, pos2) not in runs.get(iid, []):
                sys.exit('%s: %s run %d-%d is not in the .hom' % (pool, iid, pos1, pos2))
            if nsim >= len(members):
                sys.exit('%s: %s has NSIM %d in a pool of %d' % (pool, iid, nsim, len(members)))
        groups = {}
        for iid, _p1, _p2, nsim, grp in members:
            gnum = int(grp.rstrip('*'))
            groups.setdefault(gnum, []).append((iid, nsim, grp.endswith('*')))
        if sorted(groups) != list(range(1, len(groups) + 1)):
            sys.exit('%s: group numbers are not 1..n' % pool)
        for gnum, rows in groups.items():
            refs = [r for r in rows if r[2]]
            if len(refs) != 1:
                sys.exit('%s: group %d has %d references' % (pool, gnum, len(refs)))
        first = groups[1]
        ref_nsim = [r[1] for r in first if r[2]][0]
        if len(first) != ref_nsim + 1:
            sys.exit('%s: group 1 holds %d runs but its reference has NSIM %d' %
                     (pool, len(first), ref_nsim))
    print('%d pools verified' % len(pools))


if __name__ == '__main__':
    main()
