#!/usr/bin/env python3
"""Figures: manhattan enrichment plots, host-group heatmap, 3D structure map."""
import csv, sys, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

DOMAINS = {'PA': [(1, 209, 'endonuclease'), (210, 255, 'linker'), (256, 716, 'C-term (PB1/RNA binding)')],
           'PB1': [(1, 80, 'PA-binding'), (81, 600, 'RdRP core'), (601, 757, 'PB2-binding')],
           'PB2': [(1, 110, 'PB1-binding'), (250, 500, 'cap-binding region'), (530, 759, '627/NLS domain')]}

def load(path):
    return list(csv.DictReader(open(path), delimiter='\t'))

def manhattan(protein, path, outpng, markers):
    rows = load(path)
    if not rows: print('no rows', protein); return
    fig, ax = plt.subplots(figsize=(13, 4.2))
    pos = np.array([int(r['pos']) for r in rows])
    q = np.array([max(float(r['q']), 1e-300) for r in rows])
    y = -np.log10(q)
    fm = np.array([float(r['freq_m']) for r in rows])
    g2 = np.array([r['g2_pass'] == 'True' for r in rows])
    sizes = 4 + 300 * fm
    ax.scatter(pos[~g2], y[~g2], s=sizes[~g2], c='#7aa6c2', alpha=0.45, lw=0, label='tested sites')
    ax.scatter(pos[g2], y[g2], s=sizes[g2] + 10, c='#c22f2f', alpha=0.9, lw=0, label='G2 candidates')
    for (p0, p1, name) in DOMAINS.get(protein, []):
        ax.axvspan(p0, p1, color='k', alpha=0.03)
        ax.text((p0 + p1) / 2, ax.get_ylim()[1] * 0.02, name, ha='center', fontsize=7, rotation=90, color='gray')
    for k, mk in enumerate(markers):
        r = next((r for r in rows if int(r['pos']) == mk[0] and r['alt'] == mk[1]), None)
        if r:
            yv = -np.log10(max(float(r['q']), 1e-300))
            ax.annotate(mk[2], (int(r['pos']), yv), textcoords='offset points',
                        xytext=(4 + (k % 3) * 10, 5 + (k % 3) * 8), fontsize=7, color='#114',
                        arrowprops=dict(arrowstyle='-', lw=0.4, color='#889'))
    ax.axhline(-np.log10(0.05), ls='--', c='k', lw=0.7, alpha=0.6)
    ax.text(pos.max() * 0.99, -np.log10(0.05) + 0.3, 'FDR 0.05', ha='right', fontsize=7)
    ax.set_xlabel(f'{protein} position (8R1J numbering)' if protein != 'PB2' else 'PB2 position (Gs/Gd/96 numbering)')
    ax.set_ylabel('-log10(FDR q)')
    ax.set_title(f'{protein}: mammalian-enriched substitutions, H5N1 2021-2026 (point size = mammalian frequency)')
    ax.legend(loc='upper left', fontsize=7, frameon=False)
    plt.tight_layout(); plt.savefig(outpng, dpi=180); plt.close()
    print('wrote', outpng)

def heatmap(protein, group_tsv, cand_tsv, outpng):
    gf = load(group_tsv)
    if not gf: print('no groupfreq', protein); return
    groups = sorted(set(r['group'] for r in gf))
    sites = sorted(set((int(r['pos']), r['alt']) for r in gf))
    M = np.zeros((len(sites), len(groups)))
    for r in gf:
        M[sites.index((int(r['pos']), r['alt']))][groups.index(r['group'])] = float(r['freq'])
    fig, ax = plt.subplots(figsize=(7, max(2.2, 0.35 * len(sites) + 1.5)))
    im = ax.imshow(M, aspect='auto', cmap='Reds', vmin=0, vmax=max(0.2, M.max()))
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, rotation=35, ha='right', fontsize=8)
    ax.set_yticks(range(len(sites))); ax.set_yticklabels([f'{p}{a}' for p, a in sites], fontsize=8)
    for i in range(len(sites)):
        for j in range(len(groups)):
            ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=6.5,
                    color='white' if M[i, j] > 0.6 * max(0.2, M.max()) else 'black')
    ax.set_title(f'{protein} candidate frequencies by mammalian host group')
    plt.colorbar(im, label='frequency', shrink=0.8)
    plt.tight_layout(); plt.savefig(outpng, dpi=180); plt.close()
    print('wrote', outpng)

def structure3d(cand_rows, outpng):
    import sys
    sys.path.insert(0, 'code')
    from structure_map import load_coords
    c = load_coords()
    fig = plt.figure(figsize=(13, 5.6))
    for vi, (elev, azim, ttl) in enumerate([(20, -60, 'view 1'), (80, -60, 'top view')]):
        ax = fig.add_subplot(1, 2, vi + 1, projection='3d')
        for name, col, alpha in [('PB1', '#4a80b4', 0.25), ('PA', '#54a06a', 0.25), ('PB2', '#999999', 0.15), ('ANP32B', '#e08a3c', 0.8)]:
            pts = np.array(list(c[name].values()))
            ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=3, c=col, alpha=alpha, lw=0, label=name)
        for r in cand_rows:
            ch = 'PA' if r['protein'] == 'PA' else 'PB1'
            pos = int(r['pos'])
            if pos in c[ch]:
                p = c[ch][pos]
                ax.scatter([p[0]], [p[1]], [p[2]], s=60, c='#c22f2f', edgecolor='k', lw=0.4, depthshade=False)
                ax.text(p[0], p[1], p[2], f" {r['protein']}-{r['avian_cons']}{pos}{r['alt']}", fontsize=6.5)
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(ttl, fontsize=8); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        if vi == 0: ax.legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle('H5N1 polymerase + human ANP32B (PDB 8R1J, CA trace) with PB1/PA mammalian-enriched candidates', fontsize=9)
    plt.tight_layout(); plt.savefig(outpng, dpi=180); plt.close()
    print('wrote', outpng)

if __name__ == '__main__':
    mk = []
    for m in csv.DictReader(open('code/markers.tsv'), delimiter='\t'):
        mk.append((int(m['position']), m['alt_aa'], f"{m['protein']}-{m['ref_aa']}{m['position']}{m['alt_aa']}"))
    def mks(prot):
        return [(p, a, lab) for (p, a, lab) in mk if lab.startswith(prot + '-')]
    manhattan('PB1', 'results/PB1.enrichment.tsv', 'results/fig_manhattan_PB1.png', mks('PB1'))
    manhattan('PA', 'results/PA.enrichment.tsv', 'results/fig_manhattan_PA.png', mks('PA'))
    manhattan('PB2', 'results/PB2_control.enrichment.tsv', 'results/fig_manhattan_PB2.png', mks('PB2'))
    heatmap('PB1', 'results/PB1.groupfreq.tsv', None, 'results/fig_heatmap_PB1.png')
    heatmap('PA', 'results/PA.groupfreq.tsv', None, 'results/fig_heatmap_PA.png')
    cands = [r for r in load('results/PB1.enrichment.tsv') if r['g2_pass'] == 'True']
    cands += [r for r in load('results/PA.enrichment.tsv') if r['g2_pass'] == 'True']
    structure3d(cands, 'results/fig_structure3d.png')
