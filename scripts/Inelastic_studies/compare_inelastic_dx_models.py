#!/usr/bin/env python3
"""Compare simulated and polynomial inelastic dx backgrounds in GEn-II.

Independent, helicity-summed 1D template fits to the same selected events.
The polynomial models are Bernstein polynomials (degree 2, 3, ...), whose
nonnegative coefficients enforce nonnegative predicted backgrounds.

Uses uproot to stream ROOT trees, scipy for Poisson-deviance fits, and
matplotlib for comparisons/residual panels.

Run from Final_analysis repository root; see --help for options.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

DEFAULT_DATA = "/media/vimukthi/Sata/Research_JLab/outfiles/pass3/intimehighest/QE_data_GEN3_sbs100p_nucleon_np_model2_sbstrackingon_intimehighest.root"
DEFAULT_QE = "/media/vimukthi/Sata/Research_JLab/outfiles/SIMC/noNRad_correct_raster/QE_sim_simc_He3eN_GEN3_sbs100p_nucleon_np_model2_elastic_simc.root"
DEFAULT_INEL = "/media/vimukthi/Sata/Research_JLab/outfiles/inelastic_simulation_2025/QE_sim_GEN3_sbs100p_nucleon_np_model2_inelastic.root"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--data', default=DEFAULT_DATA, help='Helicity-summed experimental data ROOT file')
    parser.add_argument('--qe', default=DEFAULT_QE, help='SIMC QE p/n ROOT file')
    parser.add_argument('--inel', default=DEFAULT_INEL, help='Simulated inelastic ROOT file')
    parser.add_argument('--cuts', default='../../cuts/GEN3.json', help='Original analysis cuts JSON')
    parser.add_argument('--run-quality', default='../../DB/run_quality.json', help='Run quality JSON')
    parser.add_argument('--ignore-run-quality', action='store_true', help='Disable helicity and Moller run quality flags')
    parser.add_argument('--data-tree', default=None, help='Tree override (autodetect Tout or T)')
    parser.add_argument('--qe-tree', default=None)
    parser.add_argument('--inel-tree', default=None)
    parser.add_argument('--bins', type=int, default=200)
    parser.add_argument('--hist-range', nargs=2, type=float, default=[-2.2, 0.8], metavar=('XMIN', 'XMAX'))
    parser.add_argument('--fit-range', nargs=2, type=float, default=None, metavar=('XMIN', 'XMAX'),
                        help='Optional narrower dx fit window, same for all backgrounds')
    parser.add_argument('--degrees', nargs='+', type=int, default=[2, 3, 4, 5], help='Bernstein polynomial degrees')
    parser.add_argument('--dy-padding', type=float, default=0.05, help='dy window padding for data (and simulations in matched mode)')
    parser.add_argument('--selection-mode', choices=['current', 'matched'], default='current',
                        help='current: reproduce existing 1D data padded/simulation tight dy cuts; matched: same padded dy cut for all')
    parser.add_argument('--qe-proton-dx', type=float, default=0.0, help='Pre-shift proton simulation in meters')
    parser.add_argument('--qe-neutron-dx', type=float, default=0.0, help='Pre-shift neutron simulation in meters')
    parser.add_argument('--inel-dx', type=float, default=0.0, help='Pre-shift inelastic simulation in meters')
    parser.add_argument('--shift-range', type=float, default=0.10,
                        help='Fitted dx shift bounds +/- meters, for signal templates and sim background')
    parser.add_argument('--fix-shifts', action='store_true', help='Keep dx template shifts fixed to zero')
    parser.add_argument('--max-events', type=int, default=None,
                        help='Optional event cap PER input for quick debugging only')
    parser.add_argument('--step-size', default='100 MB', help='uproot streaming chunk size')
    parser.add_argument('--outdir', default='out_plots/inelastic_dx_background_models_GEN3', help='Output directory')
    return parser.parse_args()


def read_json(path):
    with open(path, encoding='utf-8') as stream:
        return json.load(stream)


def choose_tree(root_file, override):
    if override:
        if override not in root_file:
            raise KeyError(f"Tree '{override}' not found in {root_file.file_path}; keys={root_file.keys()}")
        return override
    for name in ('Tout', 'T'):
        if name in root_file and hasattr(root_file[name], 'iterate'):
            return name
    options = [name.split(';')[0] for name, obj in root_file.items() if hasattr(obj, 'iterate')]
    if len(options) == 1:
        return options[0]
    raise RuntimeError(f"Could not uniquely identify tree. Found: {options}. Pass --data-tree/--qe-tree/--inel-tree.")


def expressions(kind):
    common = ['dx', 'dy', 'W2', 'vz', 'eHCAL', 'ePS', 'eSH', 'trP']
    if kind == 'data':
        return common + ['runnum', 'ntrack', 'coin_time', 'helicity']
    return common + ['fnucl', 'weight'] if kind == 'qe' else common + ['weight']


def load_hist(path, tree_override, kind, cuts, quality_runs, args):
    import uproot

    bins = np.linspace(*args.hist_range, args.bins + 1)
    result = {'data': np.zeros(args.bins), 'proton': np.zeros(args.bins),
              'neutron': np.zeros(args.bins), 'inelastic': np.zeros(args.bins)}
    passed, seen = 0, 0

    with uproot.open(path) as source:
        tree_name = choose_tree(source, tree_override)
        tree = source[tree_name]
        needed = expressions(kind)
        missing = [x for x in needed if x not in tree.keys()]
        if missing:
            raise KeyError(f"{kind} tree '{tree_name}' in {path} missing {missing}. "
                           f"Available branches: {list(tree.keys())[:55]}")
        print(f'[input] {kind:>5s} {path} | tree {tree_name}, entries={tree.num_entries:,}', flush=True)
        for arrays in tree.iterate(expressions=needed, library='np', step_size=args.step_size):
            # uproot v4 and v5 normally produce a dict for flat TTree branches.
            if not isinstance(arrays, dict):
                arrays = {name: np.asarray(arrays[name]) for name in needed}
            count = len(arrays['dx'])
            if args.max_events is not None:
                available = args.max_events - seen
                if available <= 0:
                    break
                if count > available:
                    arrays = {key: value[:available] for key, value in arrays.items()}
                    count = available
            seen += count
            with np.errstate(divide='ignore', invalid='ignore'):
                ep = (arrays['ePS'] + arrays['eSH']) / arrays['trP']
            dy_pad = args.dy_padding if (kind == 'data' or args.selection_mode == 'matched') else 0.0
            selection = (np.isfinite(arrays['dx']) & np.isfinite(arrays['dy']) &
                         np.isfinite(arrays['W2']) & np.isfinite(ep) &
                         (np.abs(arrays['vz']) <= 0.27) &
                         (arrays['eHCAL'] >= cuts['eHCAL_L']) &
                         (arrays['ePS'] >= 0.2) &
                         (ep >= cuts['EoverP_L']) & (ep <= cuts['EoverP_H']) &
                         (arrays['W2'] > cuts['W2_L']) & (arrays['W2'] < cuts['W2_H']) &
                         (arrays['dy'] > cuts['dy_L'] - dy_pad) &
                         (arrays['dy'] < cuts['dy_H'] + dy_pad))
            if kind == 'data':
                selection &= ((arrays['runnum'] >= cuts['runnum_L']) &
                              (arrays['runnum'] <= cuts['runnum_H']) &
                              (arrays['ntrack'] >= 1) &
                              (np.abs(arrays['helicity']) == 1) &
                              (arrays['coin_time'] >= cuts['coin_L']) &
                              (arrays['coin_time'] <= cuts['coin_H']))
                if quality_runs is not None:
                    selection &= np.isin(arrays['runnum'], quality_runs)
                result['data'] += np.histogram(arrays['dx'][selection], bins=bins)[0]
            else:
                selection &= np.isfinite(arrays['weight'])
                if kind == 'qe':
                    p_mask = selection & (arrays['fnucl'] == 1)
                    n_mask = selection & (arrays['fnucl'] == 0)
                    result['proton'] += np.histogram(
                        arrays['dx'][p_mask] + args.qe_proton_dx,
                        bins=bins, weights=arrays['weight'][p_mask])[0]
                    result['neutron'] += np.histogram(
                        arrays['dx'][n_mask] + args.qe_neutron_dx,
                        bins=bins, weights=arrays['weight'][n_mask])[0]
                else:
                    result['inelastic'] += np.histogram(
                        arrays['dx'][selection] + args.inel_dx,
                        bins=bins, weights=arrays['weight'][selection])[0]
            passed += int(np.count_nonzero(selection))
            if args.max_events is not None and seen >= args.max_events:
                break
    print(f'[input] {kind:>5s} passed primary selection: {passed:,} of {seen:,}', flush=True)
    return result


def validate_templates(templates):
    for name, hist in templates.items():
        area = hist.sum()
        if not np.isfinite(area) or area <= 0:
            raise ValueError(f'Empty / invalid histogram: {name}; sum={area}. Check cuts and tree branches.')
        if np.any(hist < -1e-10 * abs(area)):
            raise ValueError(f'{name} histogram has substantial negative weighted bins. '
                             'Nonnegative template model would be invalid.')


def shifted_unit_shape(x, content, dx):
    """Shift bin-integrated counts using the histogram CDF (no edge spikes).

    Treat input content as constant density within each original bin. This
    gives a continuous fractional-bin translation and exactly reproduces
    the original histogram at dx=0. Normalize within the chosen fit range.
    """
    spacing = float(x[1] - x[0])
    edges = np.r_[x - 0.5 * spacing, x[-1] + 0.5 * spacing]
    cumulative = np.r_[0.0, np.cumsum(np.maximum(content, 0.0))]
    shifted_cdf = np.interp(edges - dx, edges, cumulative,
                            left=0.0, right=cumulative[-1])
    raw = np.maximum(np.diff(shifted_cdf), 0.0)
    total = raw.sum()
    return raw / total if total > 0 else np.zeros_like(raw)


def bernstein_basis(x, degree, range_lo, range_hi):
    from scipy.special import comb
    t = np.clip((x - range_lo) / (range_hi - range_lo), 0.0, 1.0)
    columns = []
    for k in range(degree + 1):
        term = comb(degree, k) * t**k * (1.0 - t)**(degree - k)
        columns.append(term / term.sum())
    return np.column_stack(columns)


def signed_poisson_deviance_residual(n, mu):
    mu = np.clip(mu, 1e-10, None)
    with np.errstate(divide='ignore', invalid='ignore'):
        contribution = mu - n + np.where(n > 0, n * np.log(np.maximum(n, 1e-300) / mu), 0.0)
    return np.sign(n - mu) * np.sqrt(2.0 * np.maximum(contribution, 0.0))


def fit_background(name, degree, x, data, p_template, n_template, i_template,
                   hist_range, shift_limit, fit_shifts):
    from scipy.optimize import least_squares
    from scipy.stats import chi2 as chi2_dist

    n_total = float(data.sum())
    shape_p = p_template / p_template.sum()
    shape_n = n_template / n_template.sum()
    shape_i = i_template / i_template.sum()
    is_sim = name == 'Simulation'
    basis = None if is_sim else bernstein_basis(x, degree, *hist_range)
    n_bg = 1 if is_sim else degree + 1
    n_shift = (3 if is_sim else 2) if fit_shifts else 0
    n_par = 2 + n_bg + n_shift
    initial = np.r_[0.65, 0.12, np.full(n_bg, 0.23 / n_bg), np.zeros(n_shift)]
    lower = np.r_[np.zeros(2 + n_bg), np.full(n_shift, -shift_limit)]
    upper = np.r_[np.full(2 + n_bg, 3.0), np.full(n_shift, shift_limit)]
    # With zero permitted shift range the fitted-shift mode would be degenerate.
    if not fit_shifts:
        assert n_shift == 0

    def components(p):
        if fit_shifts:
            dx_p, dx_n = p[2 + n_bg: 4 + n_bg]
            dx_i = p[4 + n_bg] if is_sim else 0.0
        else:
            dx_p = dx_n = dx_i = 0.0
        proton = n_total * p[0] * shifted_unit_shape(x, shape_p, dx_p)
        neutron = n_total * p[1] * shifted_unit_shape(x, shape_n, dx_n)
        if is_sim:
            background = n_total * p[2] * shifted_unit_shape(x, shape_i, dx_i)
        else:
            background = n_total * (basis @ p[2: 2 + n_bg])
        return proton, neutron, background

    def residual(parameters):
        p, n, bg = components(parameters)
        return signed_poisson_deviance_residual(data, p + n + bg)

    # Start from a few relative shifts to reduce sensitivity to local minima.
    seeds = [initial]
    if fit_shifts:
        for delta in (-0.025, 0.025):
            seed = initial.copy()
            seed[2 + n_bg] = np.clip(delta, -shift_limit, shift_limit)
            seed[3 + n_bg] = np.clip(-delta, -shift_limit, shift_limit)
            seeds.append(seed)
    outcomes = []
    for seed in seeds:
        result = least_squares(residual, seed, bounds=(lower, upper),
                               method='trf', x_scale='jac', max_nfev=2500,
                               ftol=2e-9, xtol=2e-9, gtol=2e-9)
        outcomes.append(result)
    converged = [res for res in outcomes if res.success]
    best = min(converged or outcomes, key=lambda res: np.dot(res.fun, res.fun))
    proton, neutron, background = components(best.x)
    model = proton + neutron + background
    pearson_chi2 = float(np.sum((data - model)**2 / np.clip(model, 1e-10, None)))
    deviance = float(np.sum(signed_poisson_deviance_residual(data, model)**2))
    ndf = int(data.size - n_par)
    if ndf <= 0:
        raise RuntimeError(f'Too few bins ({data.size}) for {n_par} fit parameters')
    pars = {'yield_p_frac': float(best.x[0]), 'yield_n_frac': float(best.x[1]),
            'yield_bg_frac': float(np.sum(best.x[2:2+n_bg]))}
    if is_sim:
        pars['background_coefficients'] = [float(best.x[2])]
    else:
        pars['background_coefficients'] = [float(z) for z in best.x[2:2+n_bg]]
    if fit_shifts:
        pars.update({'shift_p_m': float(best.x[2+n_bg]),
                     'shift_n_m': float(best.x[3+n_bg])})
        if is_sim:
            pars['shift_inel_m'] = float(best.x[4+n_bg])
    return dict(name=name, converged=bool(best.success), message=str(best.message),
                chi2=pearson_chi2, chi2_ndf=pearson_chi2/ndf,
                ndf=ndf, n_parameters=n_par, deviance=deviance,
                p_value=float(chi2_dist.sf(pearson_chi2, ndf)),
                aic=deviance + 2*n_par, parameters=pars,
                proton=proton, neutron=neutron, background=background,
                model=model, pulls=(data-model)/np.sqrt(np.clip(model,1e-10,None)))


def figures(x, edges, data, models, cuts, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # One consistent color scheme, similar to the ROOT components in the analysis.
    colors = {'proton': '#bd28a9', 'neutron': '#3865d6',
              'background': '#04a9aa', 'model': '#167b27'}
    names = list(models)
    dx_lo, dx_hi = cuts['dx_L'], cuts['dx_H']
    for name, fit in models.items():
        fig, (ax, pull_ax) = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                                        gridspec_kw={'height_ratios': [3.5, 1], 'hspace': 0.08})
        ax.errorbar(x, data, yerr=np.sqrt(data), fmt='.', color='black',
                    capsize=0, markersize=3, alpha=0.85, label='Data')
        for part in ('proton', 'neutron', 'background'):
            label = {'proton': 'QE proton (SIMC)', 'neutron': 'QE neutron (SIMC)',
                     'background': 'Inelastic model'}[part]
            ax.stairs(fit[part], edges, color=colors[part], lw=1.6, label=label)
        ax.stairs(fit['model'], edges, color=colors['model'], lw=1.8, label='Total fit')
        ax.axvspan(dx_lo, dx_hi, color='orange', alpha=0.09, label='Neutron dx window')
        ax.set_ylabel('Counts / bin')
        ax.set_title(f'GEN3 inelastic dx model: {name}    '
                     f'$\\chi^2_\\mathrm{{Pearson}}/\\mathrm{{ndf}}={fit["chi2"]:.1f}/{fit["ndf"]}={fit["chi2_ndf"]:.3f}$')
        ax.legend(loc='upper right', fontsize=8, ncol=2)
        ax.grid(alpha=0.12)
        pull_ax.axhline(0, lw=1.2, color='black')
        pull_ax.axhline(3, lw=0.8, color='gray', ls='--', alpha=0.6)
        pull_ax.axhline(-3, lw=0.8, color='gray', ls='--', alpha=0.6)
        pull_ax.errorbar(x, fit['pulls'], fmt='.', color='#444444', markersize=3)
        pull_ax.set_ylabel('Pull')
        pull_ax.set_xlabel(r'$dx$ (m)')
        pull_ax.grid(alpha=0.15)
        fig.subplots_adjust(left=0.09, right=0.96, bottom=0.10, top=0.91)
        filename = 'fit_simulation' if name == 'Simulation' else f'fit_{name.lower()}'
        fig.savefig(out_dir / f'{filename}.pdf')
        fig.savefig(out_dir / f'{filename}.png', dpi=180)
        plt.close(fig)

    fig, (ax, pull_ax) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                    gridspec_kw={'height_ratios':[3.2,1], 'hspace':0.08})
    ax.errorbar(x, data, yerr=np.sqrt(data), fmt='.', color='black',
                markersize=3, alpha=0.8, label='Data')
    palette = plt.get_cmap('tab10')
    for index, (name, fit) in enumerate(models.items()):
        color = palette(index)
        ax.stairs(fit['model'], edges, color=color, lw=1.7,
                  label=f'{name}: $\\chi^2_\\nu$={fit["chi2_ndf"]:.2f}')
        pull_ax.plot(x, fit['pulls'], lw=1.15, color=color, alpha=0.9, label=name)
    ax.axvspan(dx_lo, dx_hi, color='orange', alpha=0.09)
    ax.set_ylabel('Counts / bin')
    ax.set_title('Comparison of fitted inelastic-background models (same events and fit range)')
    ax.legend(fontsize=9, ncol=2)
    ax.grid(alpha=0.15)
    pull_ax.axhline(0, color='black', lw=1.1)
    pull_ax.set_ylabel('Pull')
    pull_ax.set_xlabel(r'$dx$ (m)')
    pull_ax.grid(alpha=0.15)
    fig.subplots_adjust(left=0.09, right=0.96, bottom=0.10, top=0.91)
    fig.savefig(out_dir / 'all_models_overlay.pdf')
    fig.savefig(out_dir / 'all_models_overlay.png', dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(names, [models[name]['chi2_ndf'] for name in names],
           color=[palette(i) for i in range(len(names))])
    ax.axhline(1, lw=1.2, color='black', ls='--', label=r'$\chi^2_\nu=1$')
    for i, name in enumerate(names):
        y = models[name]['chi2_ndf']
        ax.annotate(f'{y:.3f}', (i, y), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=9)
    ax.set_ylabel(r'Pearson $\chi^2/\mathrm{ndf}$')
    ax.set_title('Reduced chi-square by background model')
    ax.legend()
    ax.grid(axis='y', alpha=0.16)
    fig.tight_layout()
    fig.savefig(out_dir / 'reduced_chi2_comparison.pdf')
    fig.savefig(out_dir / 'reduced_chi2_comparison.png', dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    for i, (name, fit) in enumerate(models.items()):
        ax.stairs(fit['background'], edges, lw=1.7, label=name,
                  color=palette(i))
    ax.axvspan(dx_lo, dx_hi, color='orange', alpha=0.08)
    ax.set_xlabel(r'$dx$ (m)')
    ax.set_ylabel('Fitted inelastic-background counts / bin')
    ax.set_title('Inelastic-background shapes under the same QE proton/neutron templates')
    ax.grid(alpha=0.12)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / 'inelastic_background_overlay.pdf')
    fig.savefig(out_dir / 'inelastic_background_overlay.png', dpi=180)
    plt.close(fig)


def main():
    args = arguments()
    if args.bins < 16 or not args.hist_range[0] < args.hist_range[1]:
        raise ValueError('Use at least 16 bins and a valid --hist-range')
    if args.shift_range <= 0 and not args.fix_shifts:
        raise ValueError('Use --fix-shifts or positive --shift-range')
    if not args.degrees or any(d < 0 or d > 12 for d in args.degrees):
        raise ValueError('Polynomial degrees must be in 0..12')
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cuts = read_json(args.cuts)
    with open(out_dir / 'run_configuration.json', 'w', encoding='utf-8') as stream:
        json.dump({'command_options': vars(args), 'cuts': cuts}, stream, indent=2)
    quality_runs = None
    if not args.ignore_run_quality:
        flags = read_json(args.run_quality)
        quality_runs = np.asarray([int(run) for run, v in flags.items()
                                   if v.get('helicity', False) and v.get('moller', False)],
                                  dtype=np.int32)
        print(f'[quality] {len(quality_runs):,} runs with both flags good', flush=True)

    loaded_data = load_hist(args.data, args.data_tree, 'data', cuts, quality_runs, args)
    loaded_qe = load_hist(args.qe, args.qe_tree, 'qe', cuts, quality_runs, args)
    loaded_inel = load_hist(args.inel, args.inel_tree, 'inel', cuts, quality_runs, args)
    templates = {'data': loaded_data['data'], 'proton': loaded_qe['proton'],
                 'neutron': loaded_qe['neutron'], 'inelastic': loaded_inel['inelastic']}
    validate_templates(templates)
    print('[hist] raw 1D yields: ' + ', '.join(f'{key}={value.sum():,.2f}'
                                               for key, value in templates.items()), flush=True)
    edges = np.linspace(*args.hist_range, args.bins + 1)
    centers = 0.5*(edges[1:] + edges[:-1])
    if args.fit_range:
        select = (centers >= args.fit_range[0]) & (centers <= args.fit_range[1])
        if np.count_nonzero(select) < max(args.degrees) + 8:
            raise ValueError('Fit window too narrow for requested polynomial degrees')
        # Slice after selection: same bin subset for all models, and normalize there.
        # For contiguous selection the edges below remain correct.
        edges = edges[np.where(select)[0][0]: np.where(select)[0][-1]+2]
        centers = centers[select]
        templates = {key: value[select] for key, value in templates.items()}
        validate_templates(templates)
    print(f'[fit] x in [{edges[0]:.3f},{edges[-1]:.3f}], {centers.size} bins; dy selection mode={args.selection_mode}', flush=True)
    results = {}
    for name, degree in [('Simulation', None)] + [(f'Pol{d}', d) for d in sorted(set(args.degrees))]:
        fit = fit_background(name, degree, centers, templates['data'],
                             templates['proton'], templates['neutron'],
                             templates['inelastic'], (edges[0], edges[-1]),
                             args.shift_range, not args.fix_shifts)
        results[name] = fit
        if not fit['converged']:
            print(f'[WARNING] {name}: optimizer did not converge: {fit["message"]}', file=sys.stderr)
        print(f'[fit] {name:>10s}: chi2/ndf={fit["chi2"]:.1f}/{fit["ndf"]}='
              f'{fit["chi2_ndf"]:.4f}; deviance={fit["deviance"]:.1f}; '
              f'background={fit["parameters"]["yield_bg_frac"]:.4f} of fitted range', flush=True)
    neutron = (centers >= cuts['dx_L']) & (centers <= cuts['dx_H'])
    n_neutron_data = templates['data'][neutron].sum()
    rows = []
    for name, fit in results.items():
        predicted_bg = float(fit['background'][neutron].sum())
        row = {'model': name, 'success': fit['converged'], 'chi2': fit['chi2'],
               'ndf': fit['ndf'], 'chi2_ndf': fit['chi2_ndf'], 'p_value': fit['p_value'],
               'deviance': fit['deviance'], 'aic': fit['aic'],
               'n_parameters': fit['n_parameters'],
               'background_full_fit': float(fit['background'].sum()),
               'background_in_neutron_dx': predicted_bg,
               'data_in_neutron_dx': float(n_neutron_data),
               'background_fraction_neutron_dx': predicted_bg/n_neutron_data
                   if n_neutron_data > 0 else math.nan,
               'proton_full_fit': float(fit['proton'].sum()),
               'neutron_full_fit': float(fit['neutron'].sum())}
        rows.append(row)
    with open(out_dir / 'fit_summary.csv', 'w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with open(out_dir / 'fit_parameters.json', 'w', encoding='utf-8') as stream:
        json.dump({name: dict(converged=fit['converged'], message=fit['message'],
                              parameters=fit['parameters'], chi2=fit['chi2'],
                              ndf=fit['ndf'], aic=fit['aic'])
                   for name, fit in results.items()}, stream, indent=2)
    np.savez_compressed(out_dir / 'fit_histograms.npz', edges=edges, x=centers,
                        **{f'raw_{k}': v for k, v in templates.items()},
                        **{f'{name}_{key}': fit[key] for name, fit in results.items()
                           for key in ('proton', 'neutron', 'background', 'model', 'pulls')})
    figures(centers, edges, templates['data'], results, cuts, out_dir)
    print('\nSorted by Pearson reduced chi-square (nominal degrees of freedom):')
    print(f'{"model":>12} {"chi2/ndf":>12} {"p-value":>12} {"AIC":>12} {"bg in n-window":>16} {"bg fraction":>13}')
    for row in sorted(rows, key=lambda item: item['chi2_ndf']):
        print(f'{row["model"]:>12} {row["chi2_ndf"]:12.4f} {row["p_value"]:12.3g} '
              f'{row["aic"]:12.1f} {row["background_in_neutron_dx"]:16.1f} '
              f'{row["background_fraction_neutron_dx"]:13.5f}')
    print(f'\nWrote results to: {out_dir.resolve()}')
    print('NOTE: polynomial backgrounds may absorb signal tails; inspect plots and stability, not only reduced chi-square.')


if __name__ == '__main__':
    main()
