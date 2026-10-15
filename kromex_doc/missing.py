"""Missing chemistry (as legacy docmake): reactant sets the network never
uses (missing reactions), and products the network's reactant sets could
give but do not (missing branches). The candidates have the reactants'
elements and charge (one or two species, of the reactants' phase, no
cation-anion pair); each comes with its reaction enthalpy at 298.15 K when
the tables have every species. Gas and ice are separate."""

import itertools

MAX_SPECIES = 150      # missing reactions, per phase: they grow as N^2
MAX_BRANCH = 1000      # missing branches (Semenov's 458 gas species: 16 s)


def key_of(group, sp):
    comp, charge = {}, 0
    for x in group:
        for e, n in sp[x]["composition"].items():
            comp[e] = comp.get(e, 0) + n
        charge += sp[x]["charge"]
    return tuple(sorted(comp.items())), charge


def ion_pair(group, sp):
    c = [sp[x]["charge"] for x in group]
    return any(q > 0 for q in c) and any(q < 0 for q in c)


def repulsive(group, sp):
    c = [sp[x]["charge"] for x in group]
    return len(group) == 2 and c[0] * c[1] > 0


def analyse(net, thermo, max_species=MAX_SPECIES, max_branch=MAX_BRANCH):
    """{"phases": {phase: {"missing": [...], "branches": [...]}},
    "skipped": [notes]}"""
    sp = net.species
    out = {"phases": {}, "skipped": []}
    used = {}
    for rx in net.reactions:
        reac = tuple(sorted(r["sp"] for r in rx["reactants"]))
        prod = tuple(sorted(p["sp"] for p in rx["products"]))
        used.setdefault(reac, set()).add(prod)
    phases = {}
    for name, s in sp.items():
        if s["composition"] is None:      # opaque (grains, ...): no
            continue                      # composition to match
        phases.setdefault(s["phase"], []).append(name)
    for phase, names in sorted(phases.items()):
        names.sort()
        if len(names) > max_branch:
            out["skipped"].append(
                "%s phase: %d species, above %d; not analysed"
                % (phase, len(names), max_branch))
            continue
        do_missing = len(names) <= max_species
        if not do_missing:
            out["skipped"].append(
                "%s phase: %d species, above %d: the missing reactions "
                "are not listed (they grow as N^2), the missing branches "
                "are" % (phase, len(names), max_species))
        groups = [(a,) for a in names] + [
            (a, b) for a, b in itertools.combinations_with_replacement(
                names, 2) if not repulsive((a, b), sp)]
        by_key = {}
        for g in groups:
            by_key.setdefault(key_of(g, sp), []).append(g)

        def candidates(reac):
            return [g for g in by_key.get(key_of(reac, sp), [])
                    if g != reac and not ion_pair(g, sp)]

        def enthalpy(reac, prod):
            return thermo.dH(list(reac), list(prod))

        missing = []
        for g in groups if do_missing else ():
            if g in used:
                continue
            if len(g) == 1 and sp[g[0]]["charge"] > 0:
                continue    # a lone cation: no unimolecular process
            c = candidates(g)
            if c:
                missing.append({"reactants": list(g), "products": [
                    {"species": list(p), "dH": enthalpy(g, p)}
                    for p in c]})
        branches = []
        for reac in sorted(used):
            if not reac or any(x not in sp or sp[x]["phase"] != phase
                               or sp[x]["composition"] is None
                               for x in reac):
                continue
            have = used[reac]
            cand = [p for p in candidates(reac) if p not in have]
            if cand:
                branches.append({
                    "reactants": list(reac),
                    "present": [{"species": list(p),
                                 "dH": enthalpy(reac, p)}
                                for p in sorted(have)],
                    "missing": [{"species": list(p),
                                 "dH": enthalpy(reac, p)} for p in cand]})
        for m in missing:
            m["products"].sort(key=_order)
        for b in branches:
            b["missing"].sort(key=_order)
        out["phases"][phase] = {"missing": missing, "branches": branches,
                                "species": len(names)}
    return out


def _order(c):
    """exothermic first, then endothermic, then without data"""
    return (c["dH"] is None, c["dH"] if c["dH"] is not None else 0.0)


def mark(dH):
    """legacy's marks: exothermic, endothermic below 1e4 K, above"""
    if dH is None:
        return "no data"
    if dH < 0:
        return "exothermic"
    if dH < 1e4:
        return "endothermic < 1e4 K"
    return "endothermic"
