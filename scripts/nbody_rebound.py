#!/usr/bin/env python3
"""N-body integration of radial-velocity orbital solutions with REBOUND (WHFast).

usage: nbody_rebound.py config.txt orbits.txt summary.txt

config.txt (plain text written by Agatha's nbody.R):
    Mstar 1.0              stellar mass [Msun]
    tmax 1e7               integration time [yr]
    steps_per_orbit 20     WHFast time step = P_min / steps_per_orbit
    Nout 200               number of (log-spaced) outputs
    escape_factor 10       a particle beyond escape_factor * max(a0) has escaped
    encounter_hill 1.0     stop when two bodies come within this many mutual Hill radii
    progress_file path     optional: the driver writes 'system i n fraction' there as it goes
    system_offset 0        optional: numbering offset of the systems in the progress report
    system 1
    planet m a e omega M inc     (Msun, AU, -, rad, rad, rad; one line per planet)
    system 2
    ...

orbits.txt:  system t planet a e     (one row per output, planet and system)
summary.txt: system status t_end dE  (status: stable / escape / encounter / error)
"""
import sys, math
import rebound

def read_config(path):
    cfg = {'systems': []}
    for line in open(path):
        f = line.split()
        if not f or f[0].startswith('#'):
            continue
        if f[0] == 'system':
            cfg['systems'].append([])
        elif f[0] == 'planet':
            cfg['systems'][-1].append([float(x) for x in f[1:7]])
        elif f[0] == 'progress_file':
            cfg['progress_file'] = ' '.join(f[1:])
        else:
            cfg[f[0]] = float(f[1])
    return cfg

def report(cfg, isys, nsys, frac):
    """progress of the integration for the caller: 'system i n fraction'"""
    pf = cfg.get('progress_file')
    if pf:
        try:
            with open(pf, 'w') as fh:
                fh.write('%d %d %.4f\n' % (isys, nsys, frac))
        except OSError:
            pass

def hill_min(sim):
    """the smallest mutual Hill radius of adjacent planets"""
    ps = sorted(sim.particles[1:], key=lambda p: p.a)
    M = sim.particles[0].m
    rh = []
    for p1, p2 in zip(ps[:-1], ps[1:]):
        rh.append(((p1.m + p2.m) / (3.0 * M)) ** (1.0 / 3.0) * 0.5 * (p1.a + p2.a))
    return min(rh) if rh else 0.0

def main():
    cfg = read_config(sys.argv[1])
    orb = open(sys.argv[2], 'w')
    summ = open(sys.argv[3], 'w')
    orb.write('system t planet a e\n')
    summ.write('system status t_end dE\n')
    tmax = cfg.get('tmax', 1e7)
    Nout = int(cfg.get('Nout', 200))
    off = int(cfg.get('system_offset', 0))
    ntot = int(cfg.get('system_total', len(cfg['systems'])))
    for isys, planets in enumerate(cfg['systems'], start=1):
        report(cfg, isys + off, ntot, 0.0)
        sim = rebound.Simulation()
        sim.units = ('yr', 'AU', 'Msun')
        sim.add(m=cfg.get('Mstar', 1.0))
        for m, a, e, omega, M, inc in planets:
            sim.add(m=m, a=a, e=e, omega=omega, M=M, inc=inc)
        sim.move_to_com()
        sim.integrator = 'whfast'
        Pmin = min(p.P for p in sim.particles[1:])
        sim.dt = Pmin / cfg.get('steps_per_orbit', 20)
        amax0 = max(p.a for p in sim.particles[1:])
        sim.exit_max_distance = cfg.get('escape_factor', 10.0) * amax0
        rh = hill_min(sim)
        if rh > 0 and cfg.get('encounter_hill', 1.0) > 0:
            sim.exit_min_distance = cfg.get('encounter_hill', 1.0) * rh
        E0 = sim.energy()
        t1 = max(10.0 * Pmin, tmax / 1e4)
        times = [0.0] + [t1 * (tmax / t1) ** (k / (Nout - 1.0)) for k in range(Nout)] if tmax > t1 else [0.0, tmax]
        status, t_end = 'stable', tmax
        try:
            for t in times:
                sim.integrate(t, exact_finish_time=0)
                report(cfg, isys + off, ntot, sim.t / tmax)
                bad = False
                for j, o in enumerate(sim.orbits(), start=1):
                    orb.write('%d %.9g %d %.9g %.9g\n' % (isys, sim.t, j, o.a, o.e))
                    if not (o.a > 0 and o.e < 1):
                        bad = True
                if bad:
                    status, t_end = 'escape', sim.t
                    break
        except rebound.Escape:
            status, t_end = 'escape', sim.t
        except rebound.Encounter:
            status, t_end = 'encounter', sim.t
        except Exception as err:
            status, t_end = 'error', sim.t
            sys.stderr.write('system %d: %s\n' % (isys, err))
        dE = abs((sim.energy() - E0) / E0) if E0 != 0 else float('nan')
        summ.write('%d %s %.9g %.3g\n' % (isys, status, t_end, dE))
        print('system %d/%d %s t=%.4g yr' % (isys, len(cfg['systems']), status, t_end), flush=True)
    orb.close(); summ.close()

if __name__ == '__main__':
    main()
