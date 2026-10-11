#!/usr/bin/env python3
"""Translate CCINF truncations between lattice and physical units, and size jobs.

Forward:  given Nmax / t3_cut, print kmax, kmax/kF, the 3p3h cut in MeV and in
          units of e_F, the orbital count, the T3 dimension and memory.
Inverse:  given --kmax/--kmax-kf and --e3/--e3-ef, pick the Nmax and t3_cut
          that reach them at each (A, rho), then print the same table.

The counting mirrors the Fortran: basis.f90 (orbitals, holes, NF2),
modules.f90 cut3b / number_3b_confs, and ndim3 in setup_t3_amplitudes.

Examples:
  ccinf_params.py pnm 66 0.16 --nmax 10 --t3cut 12
  ccinf_params.py snm 132 228 --rho 0.08 0.12 0.16 --kmax 3.0 --e3 300 --ranks 256
"""
import argparse
import math
import numpy as np

HBARC = 197.3269804
NUC_MASS = 0.5 * (938.27208816 + 939.56542052)
NOCCLIST = [0, 2, 14, 38, 54, 66, 114, 162, 186, 246]          # per species, index = Nocc
NMAXLIST = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16,    # index = Nmax - 1
            17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30,
            32, 33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44]


def species_tz(species):
    return [1] if species == 'pnm' else [1, -1]


def per_species_count(species, A):
    n = A if species == 'pnm' else A // 2
    if n not in NOCCLIST or (species == 'snm' and A % 2):
        ok = NOCCLIST[1:] if species == 'pnm' else [2 * x for x in NOCCLIST[1:]]
        raise SystemExit(f'A={A} is not a closed shell for {species}; use one of {ok}')
    return n


def scales(species, A, rho):
    L = (A / rho) ** (1 / 3)
    rho_s = rho if species == 'pnm' else rho / 2
    kf = (3 * math.pi**2 * rho_s) ** (1 / 3)
    eps1 = HBARC**2 / (2 * NUC_MASS) * (2 * math.pi / L) ** 2   # MeV per unit of n^2
    ef = HBARC**2 * kf**2 / (2 * NUC_MASS)
    return L, kf, eps1, ef


def pick_nmax(kmax, L):
    n2 = (kmax * L / (2 * math.pi)) ** 2
    for i, v in enumerate(NMAXLIST):
        if v >= n2 - 1e-9:
            return i + 1
    raise SystemExit(f'kmax={kmax:.3f} needs n^2={n2:.1f} > {NMAXLIST[-1]} (end of Nmaxlist)')


def pair_tables(ga, gb, same, M, Wmax):
    """Orbital pairs a<b between group sets ga, gb (spin multiplicity folded in).

    Returns per-channel pair counts and cumulative counts by summed weight,
    with channels laid out on a (4M+1)^3 x 3(Tz) grid.
    """
    na, ta, wa = ga
    nb, tb, wb = gb
    if same:
        i, j = np.triu_indices(len(na))
        mult = np.where(i == j, 1, 4)
    else:
        i, j = np.meshgrid(np.arange(len(na)), np.arange(len(nb)), indexing='ij')
        i, j = i.ravel(), j.ravel()
        mult = np.full(i.size, 4)
    N = na[i] + nb[j] + 2 * M
    Tz = (ta[i] + tb[j]) // 2 + 1
    W = wa[i] + wb[j]
    side = 4 * M + 1
    ch = ((N[:, 0] * side + N[:, 1]) * side + N[:, 2]) * 3 + Tz
    nch = side**3 * 3
    count = np.bincount(ch, weights=mult, minlength=nch)
    hist = np.zeros((nch, Wmax + 1))
    np.add.at(hist, (ch, W), mult)
    cum = np.cumsum(hist, axis=1).reshape(side, side, side, 3, Wmax + 1)
    return count, cum


def three_body_counts(cum, groups, cut, M, Wmax):
    """number_3b_confs for every 3b channel: sum over c of pairs in channel (N3-n_c)
    whose weight leaves room for w_c under cut3b."""
    side3 = 6 * M + 1
    out = np.zeros((side3, side3, side3, 4))   # tz3 = -3,-1,1,3
    n, tz, w = groups
    for nc, tc, wc in zip(n, tz, w):
        thr = cut + 3 - wc
        if thr < 0:
            continue
        thr = min(thr, Wmax)
        sl = cum[..., thr]                      # (4M+1)^3 x Tz1
        ox, oy, oz = nc + M
        for t1 in range(3):
            tz3 = 2 * (t1 - 1) + tc
            out[ox:ox + 4 * M + 1, oy:oy + 4 * M + 1, oz:oz + 4 * M + 1, (tz3 + 3) // 2] += 2 * sl[..., t1]
    return out


# Riker cost model, fitted to the 2026-10-10 trial (nnlosat, rho=0.16, tnf=1, 16 ranks x 8 threads/node):
# CCD Nmax 10-30 for pnm66/snm132 on 4 nodes, CCDT-1 snm132 Nmax 16/20 (4 nodes) and 24 (8 nodes).
NODE_GB = 2200.0
RANKS_PER_NODE = 16
MEM_REP = (56.8, 0.083)          # GB/rank = a*t2_gb + b*pppp_gb (replicated, rank-independent; +-8%)
T_CCD = (5.89e-7, 2.0, 1.31)     # CCD seconds on 4 nodes = c * orbs^p * holes^q (+-8% rms)
T_T3 = (200.0, 6.7e-9)           # CCDT-1 extra seconds = a + b * ndim3 * 4/nodes (+-20-30%)


def cost(row, nodes, mem_frac):
    """Memory per rank, minimum nodes and wall-time estimates for one analysed row."""
    rep = MEM_REP[0] * row['pphh_gb'] + MEM_REP[1] * row['pppp_gb']
    t3 = row['t3_gb'] * 1.01                      # + antisym slab / imbalance (code's peak is ~0.3% over)
    room = mem_frac * NODE_GB - RANKS_PER_NODE * rep
    nmin = math.inf if room <= 0 else max(1, math.ceil(t3 / room))
    n = nodes if nodes else nmin
    gb_rank = rep + t3 / (RANKS_PER_NODE * n) if n != math.inf else math.inf
    t_ccd = T_CCD[0] * row['orbs']**T_CCD[1] * row['holes']**T_CCD[2]
    t_t3 = T_T3[0] + T_T3[1] * row['ndim3'] * 4 / n if n != math.inf else math.inf
    return dict(rep_gb=rep, nodes_min=nmin, nodes=n, gb_rank=gb_rank,
                node_frac=RANKS_PER_NODE * gb_rank / NODE_GB, t_ccd=t_ccd, t_ccdt1=t_ccd + t_t3)


NMAXLIST0 = [0, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14,   # index = n^2 (get_eom_ind)
             15, 16, 17, 18, 19, 20, 20, 21, 22, 23, 24, 24, 25, 26,
             26, 27, 28, 29, 30, 30, 32, 33, 34, 35, 36, 37, 38, 39]
EOM_ITERATIONS = 550


def code_order(n2max, tzs):
    """Orbitals in setup_sp_data's loop order: shell, max|n_i|, |nx|,|ny|,|nz|, signs, sz, tz."""
    M = math.isqrt(n2max)
    mags = [(x, y, z) for x in range(M + 1) for y in range(M + 1) for z in range(M + 1)
            if x * x + y * y + z * z <= n2max]
    mags.sort(key=lambda m: (m[0]**2 + m[1]**2 + m[2]**2, max(m), m[0], m[1], m[2]))
    out = []
    for x, y, z in mags:
        for ix in ([-1, 1] if x else [0]):
            for iy in ([-1, 1] if y else [0]):
                for iz in ([-1, 1] if z else [0]):
                    for _sz in (-1, 1):
                        for tz in sorted(tzs):
                            out.append(((ix * x, iy * y, iz * z), tz))
    return out


def eom_targets(n2max, tzs, nf2, tz0=1):
    """get_eom_ind for add = -4..-1, 1..4 (neutron, as with add_n=100)."""
    order = code_order(n2max, tzs)
    mine = [(n, sum(c * c for c in n)) for n, tz in order if tz == tz0]
    holes = [m for m in mine if m[1] <= nf2]
    parts = [m for m in mine if m[1] > nf2]
    targets = {}
    for add in (1, 2, 3, 4):
        targets[add] = next((n for n, s in parts if NMAXLIST0[s] - NMAXLIST0[nf2] >= add), None)
    low = parts[0][1] if parts else None
    for add in (-1, -2, -3, -4):
        targets[add] = None if low is None else \
            next((n for n, s in reversed(holes) if NMAXLIST0[s] - NMAXLIST0[low] <= add), None)
    return targets


def degenerate_counts(groups, M):
    """Ordered pairs of distinct orbitals (x, y) on the 3b grid with N3 = 2n_x + n_y:
    the (c, ab) entries of number_3b_confs where c coincides with a or b."""
    n, tz, _ = groups
    g, h = np.meshgrid(np.arange(len(n)), np.arange(len(n)), indexing='ij')
    g, h = g.ravel(), h.ravel()
    mult = np.where(g == h, 2, 4)
    N = 2 * n[g] + n[h] + 3 * M
    t = (2 * tz[g] + tz[h] + 3) // 2
    side3 = 6 * M + 1
    idx = ((N[:, 0] * side3 + N[:, 1]) * side3 + N[:, 2]) * 4 + t
    return np.bincount(idx, weights=mult, minlength=side3**3 * 4).reshape(side3, side3, side3, 4)


def eom_sizes(P, H, pp, hh, M, Wmax, cum_pp, cum_hh, targets, eom3cut):
    """EOM2/EOM3 dimensions per target, as setup_pa/preom(3)_amplitudes allocate them.

    The EOM3 amplitudes and Krylov vectors run over all pairs of each channel: eom3_cut
    only selects which 3b channels exist (mappings.f90), so 'alloc' and 'vec' are uncut.
    """
    side = 4 * M + 1
    pp = pp.reshape(side, side, side, 3)
    hh = hh.reshape(side, side, side, 3)
    big = 10**6
    ppp_cut = three_body_counts(cum_pp, P, eom3cut, M, Wmax)
    hhh_cut = three_body_counts(cum_hh, H, eom3cut, M, Wmax)
    ppp_all = three_body_counts(cum_pp, P, big, M, Wmax)
    hhh_all = three_body_counts(cum_hh, H, big, M, Wmax)
    ppp_tri = (ppp_all - degenerate_counts(P, M)) / 3
    hhh_tri = (hhh_all - degenerate_counts(H, M)) / 3

    def two_body_at(grid, n3off, tzoff):
        """grid[N3 - n_t, (tz3 - tz_t)/2] laid out on the 3b grid (zero where out of range)."""
        side3 = 6 * M + 1
        out = np.zeros((side3, side3, side3, 4))
        for t3i in range(4):
            tz2 = (2 * t3i - 3 - tzoff)
            if tz2 % 2 or abs(tz2 // 2) > 1:
                continue
            o = [n3off[k] + M for k in range(3)]
            out[o[0]:o[0] + side, o[1]:o[1] + side, o[2]:o[2] + side, t3i] = grid[..., tz2 // 2 + 1]
        return out

    def scanned(groups):
        """setup_pa/preom3 only scan N3 in [n3min, n3max] (every component) and tz3 in
        [t3min, t3max], taken from triples of three distinct orbitals in the set."""
        nx = np.sort(np.repeat(groups[0][:, 0], 2))
        tz = np.sort(np.repeat(groups[1], 2))
        side3 = 6 * M + 1
        mask = np.zeros((side3, side3, side3, 4), dtype=bool)
        if len(nx) < 3:
            return mask
        lo, hi = nx[:3].sum() + 3 * M, nx[-3:].sum() + 3 * M
        tlo, thi = (tz[:3].sum() + 3) // 2, (tz[-3:].sum() + 3) // 2
        mask[lo:hi + 1, lo:hi + 1, lo:hi + 1, tlo:thi + 1] = True
        return mask

    in_pa, in_pr = scanned(P), scanned(H)

    def pairs_with(grid, nt, groups):
        """sum over orbitals x in groups of grid[nt + n_x, (1 + tz_x)/2] (neutron target)."""
        idx = groups[0] + nt + 2 * M
        return int(2 * grid[idx[:, 0], idx[:, 1], idx[:, 2], (1 + groups[1]) // 2 + 1].sum())

    res = {}
    for add, nt in targets.items():
        if nt is None:
            continue
        nt = np.array(nt)
        if add > 0:     # PA: 1p + 2p1h, then 3p2h
            n2b = 2 + pairs_with(pp, nt, H)
            other = two_body_at(hh, nt, 1)
            sel = (ppp_cut > 0) & (other > 0) & in_pa
            alloc = float((ppp_all * other)[sel].sum())
            vec = float((ppp_tri * other)[sel].sum())
            prt = float((ppp_cut * other)[sel].sum())
        else:           # PR: 1h + 2h1p, then 3h2p
            n2b = 2 + pairs_with(hh, nt, P)
            other = two_body_at(pp, nt, 1)
            sel = (hhh_cut > 0) & (other > 0) & in_pr
            alloc = float((hhh_all * other)[sel].sum())
            vec = float((hhh_tri * other)[sel].sum())
            prt = float((hhh_cut * other)[sel].sum())
        res[add] = dict(target=tuple(int(x) for x in nt), ndim2=int(n2b), printed3=int(prt),
                        alloc3=int(alloc), vec3=int(round(vec)),
                        eom2_gb=(EOM_ITERATIONS + 2) * n2b * 16 / 1e9,
                        eom3_gb=((EOM_ITERATIONS + 2) * (n2b + vec) + alloc) * 16 / 1e9)
    return res


def analyse(species, A, rho, nmax, t3cut, tnf, ranks, count, eom3cut=None):
    nper = per_species_count(species, A)
    L, kf, eps1, ef = scales(species, A, rho)
    n2max = NMAXLIST[nmax - 1]
    M = math.isqrt(n2max)
    r = np.arange(-M, M + 1)
    vec = np.array(np.meshgrid(r, r, r, indexing='ij')).reshape(3, -1).T
    vec = vec[(vec**2).sum(1) <= n2max]
    n2 = (vec**2).sum(1)

    # holes fill whole shells: nF^2 is the shell where the per-species count reaches nper
    shells = np.unique(n2)
    filled = np.cumsum([2 * (n2 == s).sum() for s in shells])
    if nper not in filled:
        raise SystemExit(f'Nmax={nmax} (n^2<={n2max}) cannot hold {nper} particles per species')
    nf2 = shells[list(filled).index(nper)]

    tzs = species_tz(species)
    n = np.concatenate([vec] * len(tzs))
    tz = np.repeat(tzs, len(vec))
    d = np.concatenate([n2] * len(tzs)) - nf2
    w = np.where(d <= 0, np.abs(d) + 1, d)
    hole = d <= 0
    orbs = 2 * len(n)

    row = dict(species=species, A=A, rho=rho, L=L, kf=kf, nmax=nmax, n2max=n2max,
               kmax=2 * math.pi * math.sqrt(n2max) / L, t3cut=t3cut,
               e3=(t3cut + 3) * eps1, ef=ef, orbs=orbs, holes=2 * hole.sum())
    if not count:
        return row

    P = (n[~hole], tz[~hole], w[~hole])
    H = (n[hole], tz[hole], w[hole])
    Wmax = 2 * int(w.max())
    pp, cum_pp = pair_tables(P, P, True, M, Wmax)
    hh, cum_hh = pair_tables(H, H, True, M, Wmax)
    ppp = three_body_counts(cum_pp, P, t3cut, M, Wmax)
    hhh = three_body_counts(cum_hh, H, t3cut, M, Wmax)
    ndim3 = int((ppp * hhh).sum())
    t3_bytes = 8                                  # complex(t3c), t3c = single precision
    nstore = 2 if tnf > 1 else 1
    row.update(ndim3=ndim3,
               t3_gb=nstore * ndim3 * t3_bytes / 1e9,
               pppp_gb=float((pp**2).sum()) * 16 / 1e9,
               pphh_gb=float((pp * hh).sum()) * 16 / 1e9)
    row['t3_rank_gb'] = row['t3_gb'] / ranks
    if eom3cut is not None:
        targets = eom_targets(n2max, tzs, nf2)
        row['eom'] = eom_sizes(P, H, pp, hh, M, Wmax, cum_pp, cum_hh, targets, eom3cut)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('species', choices=['pnm', 'snm'])
    ap.add_argument('A', type=int, nargs='+', help='particle number(s), closed shells')
    ap.add_argument('--rho', type=float, nargs='+', default=[0.16], help='density/densities in fm^-3')
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--nmax', type=int, help='Nmax input of ccdt.exe')
    g.add_argument('--kmax', type=float, help='target momentum cutoff in fm^-1')
    g.add_argument('--kmax-kf', type=float, help='target momentum cutoff in units of kF')
    g3 = ap.add_mutually_exclusive_group()
    g3.add_argument('--t3cut', type=int, help='t3_cut input of ccdt.exe (default: no cut)')
    g3.add_argument('--e3', type=float, help='target 3p3h cut in MeV')
    g3.add_argument('--e3-ef', type=float, help='target 3p3h cut in units of e_F')
    ap.add_argument('--tnf', type=int, default=2, help='tnf_approx; 2 also stores t3_ccm0 (default 2)')
    ap.add_argument('--ranks', type=int, default=1, help='MPI ranks for the per-rank T3 column')
    ap.add_argument('--no-count', action='store_true', help='skip the T3/2b counting (fast)')
    ap.add_argument('--eom3cut', type=int, help='also size the add_n=100 EOM states with this eom3_cut '
                    '(1000 = none); EOM3 is freed T3, so the run peak is max(T3, EOM3)')
    ap.add_argument('--cost', action='store_true', help='Riker estimate: GB/rank, minimum nodes, CCD and '
                    'CCDT-1 wall time at 16 ranks x 8 threads per node')
    ap.add_argument('--nodes', type=int, help='with --cost: use this node count instead of the minimum')
    ap.add_argument('--mem-frac', type=float, default=0.8, help='with --cost: usable fraction of node memory '
                    '(default 0.8)')
    a = ap.parse_args()
    if a.cost and a.no_count:
        ap.error('--cost needs the T3/2b counting; drop --no-count')

    hdr = (f"{'sp':>3} {'A':>4} {'rho':>5} {'L':>6} {'kF':>5} | {'Nmax':>4} {'n2':>3} {'kmax':>5} {'k/kF':>5} | "
           f"{'t3cut':>5} {'E3':>7} {'E3/eF':>6} | {'orbs':>5} {'holes':>5}")
    if not a.no_count:
        hdr += f" | {'T3 dim':>14} {'T3 GB':>9} {'GB/rank':>8} {'pppp GB':>8} {'t2 GB':>7}"
    print(hdr)
    print('-' * len(hdr))
    for A in a.A:
        for rho in a.rho:
            L, kf, eps1, ef = scales(a.species, A, rho)
            if a.nmax is not None:
                nmax = a.nmax
            else:
                nmax = pick_nmax(a.kmax if a.kmax is not None else a.kmax_kf * kf, L)
            if a.t3cut is not None:
                t3cut = a.t3cut
            elif a.e3 is not None or a.e3_ef is not None:
                e3 = a.e3 if a.e3 is not None else a.e3_ef * ef
                t3cut = max(0, math.floor(e3 / eps1 + 1e-9) - 3)
            else:
                t3cut = 1000
            r = analyse(a.species, A, rho, nmax, t3cut, a.tnf, a.ranks, not a.no_count,
                        None if a.no_count else a.eom3cut)
            cut = 'none' if t3cut >= 100 else f'{t3cut:5d}'
            e3s = '-' if t3cut >= 100 else f"{r['e3']:7.1f}"
            e3f = '-' if t3cut >= 100 else f"{r['e3'] / r['ef']:6.2f}"
            line = (f"{r['species']:>3} {A:4d} {rho:5.3f} {r['L']:6.2f} {r['kf']:5.3f} | {nmax:4d} {r['n2max']:3d} "
                    f"{r['kmax']:5.2f} {r['kmax'] / r['kf']:5.2f} | {cut:>5} {e3s:>7} {e3f:>6} | "
                    f"{r['orbs']:5d} {r['holes']:5d}")
            if not a.no_count:
                line += (f" | {r['ndim3']:14d} {r['t3_gb']:9.2f} {r['t3_rank_gb']:8.3f} "
                         f"{r['pppp_gb']:8.2f} {r['pphh_gb']:7.3f}")
            print(line)
            if a.cost:
                c = cost(r, a.nodes, a.mem_frac)
                if c['nodes'] == math.inf:
                    print(f"      cost: replicated {c['rep_gb']:.1f} GB/rank x {RANKS_PER_NODE} exceeds "
                          f"{a.mem_frac:.0%} of a node")
                else:
                    warn = '  (over --mem-frac)' if c['node_frac'] > a.mem_frac else ''
                    print(f"      cost: {c['nodes']} nodes (min {c['nodes_min']}) x {RANKS_PER_NODE}x8 | "
                          f"{c['gb_rank']:.1f} GB/rank ({c['rep_gb']:.1f} replicated), {c['node_frac']:.0%} of node"
                          f"{warn} | CCD {c['t_ccd'] / 60:.0f} min, CCDT-1 {c['t_ccdt1'] / 60:.0f} min")
            if 'eom' in r:
                print(f"      {'add':>4} {'target n':>12} {'EOM2 dim':>10} {'EOM3 dim*':>14} {'EOM3 vec':>14} "
                      f"{'EOM3 alloc':>14} {'EOM2 GB':>8} {'EOM3 GB':>9}")
                for add in sorted(r['eom']):
                    e = r['eom'][add]
                    print(f"      {add:4d} {str(e['target']):>12} {e['ndim2']:10d} {e['printed3']:14d} {e['vec3']:14d} "
                          f"{e['alloc3']:14d} {e['eom2_gb']:8.3f} {e['eom3_gb']:9.2f}")
                eom3 = max(e['eom3_gb'] for e in r['eom'].values())
                print(f"      peak with EOM3: max(T3 {r['t3_gb']:.2f}, EOM3 {eom3:.2f}) = {max(r['t3_gb'], eom3):.2f} GB")
    print('\nkmax = 2pi*sqrt(n2)/L; E3 = (t3_cut+3)*hbar^2(2pi/L)^2/2m, the 3-particle excitation above the')
    print('Fermi surface allowed by cut3b. T3 GB is single precision, x2 with t3_ccm0 when tnf=2; GB/rank')
    print('assumes perfect balance and excludes the antisym temp. pppp is distributed, t2 (pphh) replicated.')
    if a.cost:
        print('cost: fit to the 2026-10-10 Riker trial (A=66 pnm, 132 snm, rho=0.16, tnf=1). Memory +-8%; CCD time')
        print('+-8% (node count assumed not to matter); CCDT-1 extra +-20-30% (iteration counts vary). Other A')
        print('extrapolate through orbs^2*holes^1.31: check one point before trusting. tnf=2 V3 fill not timed.')
    if a.eom3cut is not None:
        print('EOM: one state at a time; GB = 552 Krylov vectors (+ r3 for EOM3) in double complex, Hbar excluded.')
        print('EOM3 dim* is the cut count the code prints; the code allocates "alloc" and iterates "vec",')
        print('both over all pairs, because eom3_cut only picks which 3b channels exist (mappings.f90).')


if __name__ == '__main__':
    main()
