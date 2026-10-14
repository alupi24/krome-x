"""The generation options of KROME-X: one definition per option, for the
Python object (kromex.Preprocessor), the TOML file and --set
(docs/config-design.md).

Each option is a section and a key (thermal.cooling); the definition holds
its type, its allowed values, its default and the generator argument it
sets. Values are checked here, when set; the rules that tie options
together (thermal.metals needs Z in thermal.cooling, ...) stay in the
generator."""

import difflib


class ConfigError(ValueError):
    """an unknown option, or a value it does not accept"""


class Opt:
    """one option. kind: path, str, enum, int, float, bool, notbool (true
    sets the generator's negative flag to false), list, floor (false, true
    or a list) or bindings (a list of fortran, python)."""
    __slots__ = ("section", "key", "dest", "kind", "choices", "items",
                 "default", "flag", "open")

    def __init__(self, section, key, dest, kind, default, choices=None,
                 items=None, flag=None, open=False):
        self.section, self.key, self.dest, self.kind = section, key, dest, kind
        self.default, self.choices, self.items = default, choices, items
        # open: items outside the known ones are passed to the generator
        # (upper case), which knows them (the coolants of coolant files)
        self.open = open
        # the generator's command-line flag (the internal bridge to it)
        self.flag = flag or ("--" + dest.replace("_", "-") if dest else None)

    @property
    def name(self):
        return "%s.%s" % (self.section, self.key)


def _codegen():
    import codegen  # noqa: the generator, for its tables of names
    return codegen


def _cooling():
    """the built-in cooling processes, and FILE (every coolant of the
    coolant files); a coolant file's own names are accepted too (open)"""
    return tuple(_codegen().Codegen.COOLING_NAMES) + ("FILE",)


def _heating():
    return tuple(_codegen().Codegen.HEATING_NAMES)


def _floor():
    return tuple(_codegen().Codegen.FLOOR_NAMES)


# the sections, in the order of the template and of manual section 3
SECTIONS = ("project", "rates", "thermal", "radiation", "dust", "surface",
            "solver")

OPTIONS = (
    # project
    Opt("project", "network", "network", "path", None),
    Opt("project", "outdir", "outdir", "path", None),
    Opt("project", "bindings", None, "bindings", [],
        items=("fortran", "python")),
    Opt("project", "single", "single", "bool", False),
    Opt("project", "runtime", "no_runtime", "notbool", True),
    # rates
    Opt("rates", "outside", "outside", "enum", "asym",
        choices=("asym", "clamp", "zero")),
    Opt("rates", "rhs", "rhs", "enum", "auto",
        choices=("auto", "unrolled", "loop")),
    Opt("rates", "reactant_multiplicity", "reactant_multiplicity", "bool",
        False),
    Opt("rates", "screening", "screening", "bool", False),
    Opt("rates", "check_reverse", "check_reverse", "bool", False),
    # thermal
    Opt("thermal", "cooling", "cooling", "list", [], items=_cooling,
        open=True),
    Opt("thermal", "heating", "heating", "list", [], items=_heating),
    Opt("thermal", "floor", "cooling_floor", "floor", False, items=_floor),
    Opt("thermal", "gamma", "gamma", "enum", "CONST",
        choices=("CONST", "EXACT", "POPOVAS")),
    Opt("thermal", "h2_opacity", "h2_opacity", "enum", "",
        choices=("", "RIPAMONTI", "OMUKAI")),
    Opt("thermal", "metals", "metal_cooling", "list", []),
    Opt("thermal", "coolant_files", "cool_file", "list", []),
    Opt("thermal", "escape", "cool_escape", "enum", "",
        choices=("", "LVG")),
    Opt("thermal", "levels", "cool_levels", "int", None),
    Opt("thermal", "photo_pumping", "photo_induced", "bool", False),
    Opt("thermal", "cmb_pumping", "cmb_pumping", "bool", False),
    Opt("thermal", "photoelectric", "photoelectric", "enum", "BT94",
        choices=("BT94", "WD01")),
    Opt("thermal", "wd01", "wd01", "list", []),
    # radiation
    Opt("radiation", "evolve", "evolve_rad", "bool", False),
    Opt("radiation", "shielding", "shielding", "enum", "WG11",
        choices=("WG11", "DB96", "R14", "NONE")),
    Opt("radiation", "num2col", "num2col", "enum", "default",
        choices=("default", "jeans", "jeans40", "length")),
    Opt("radiation", "habing_dust_shielding", "habing_dust_shielding",
        "bool", False),
    Opt("radiation", "secondary_ionization", "secondary_ionization", "list",
        []),
    # dust
    Opt("dust", "tables", "dust_tables", "enum", "",
        choices=("", "HM2012", "MC", "3D", "PHOTO3D")),
    Opt("dust", "evolve_Td", "evolve_Td", "bool", False,
        flag="--evolve-Td"),
    Opt("dust", "evolve", "evolve_dust", "list", [],
        items=("GROWTH", "SPUTTER", "EVAP", "KELVIN")),
    Opt("dust", "gas_total", "gas_total", "bool", False),
    # surface
    Opt("surface", "binding_energies", "binding_energies", "enum", "KROME",
        choices=("KROME", "SEMENOV")),
    Opt("surface", "nu", "surface_nu", "enum", "HASEGAWA",
        choices=("HASEGAWA", "LEGACY")),
    Opt("surface", "competition", "no_surface_competition", "notbool", True),
    Opt("surface", "tunnelling", "no_surface_tunnelling", "notbool", True),
    Opt("surface", "reactive_desorption", "reactive_desorption", "enum",
        "NONE", choices=("NONE", "FIXED", "MINISSALE")),
    Opt("surface", "bins", "surface_bins", "int", 0),
    Opt("surface", "cr_desorption", "cr_desorption", "str", ""),
    # solver defaults (the caller can still change them at run time)
    Opt("solver", "rtol", "rtol", "float", None),
    Opt("solver", "atol", "atol", "float", None),
    Opt("solver", "atol_file", "custom_atol", "path", None),
    Opt("solver", "rtol_file", "custom_rtol", "path", None),
)

INDEX = {o.name: o for o in OPTIONS}


def section_keys(section):
    return [o.key for o in OPTIONS if o.section == section]


def lookup(name):
    """the option SECTION.KEY, or a ConfigError suggesting close names"""
    if name in INDEX:
        return INDEX[name]
    sec = name.split(".", 1)[0]
    close = difflib.get_close_matches(name, list(INDEX), n=3, cutoff=0.6)
    if sec not in SECTIONS:
        close = close or difflib.get_close_matches(sec, SECTIONS, n=2)
        raise ConfigError("unknown section in %r (sections: %s)%s"
                          % (name, ", ".join(SECTIONS), _did(close)))
    raise ConfigError("unknown option %r (the %s options: %s)%s"
                      % (name, sec, ", ".join(section_keys(sec)),
                         _did(close)))


def _did(close):
    return "; did you mean %s?" % " or ".join(close) if close else ""


def _items(opt):
    it = opt.items
    return it() if callable(it) else it


def _pick(opt, value, allowed, what):
    """value matched case-insensitively against allowed, as spelled there"""
    for a in allowed:
        if str(value).upper() == str(a).upper():
            return a
    close = difflib.get_close_matches(str(value), [str(a) for a in allowed],
                                      n=2, cutoff=0.5)
    raise ConfigError("%s: %s %r is not valid (valid: %s)%s"
                      % (opt.name, what, value,
                         ", ".join(repr(a) if a == "" else str(a)
                                   for a in allowed), _did(close)))


def check(opt, value):
    """value checked and normalized for opt (ConfigError otherwise)"""
    k = opt.kind
    if value is None:
        if opt.default is None:
            return None
        raise ConfigError("%s cannot be None" % opt.name)
    if k in ("bool", "notbool"):
        if not isinstance(value, bool):
            raise ConfigError("%s takes true or false, not %r"
                              % (opt.name, value))
        return value
    if k == "int":
        if isinstance(value, bool) or not isinstance(value, int) \
                or value < 0:
            raise ConfigError("%s takes a non-negative integer, not %r"
                              % (opt.name, value))
        return value
    if k == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)) \
                or value <= 0:
            raise ConfigError("%s takes a positive number, not %r"
                              % (opt.name, value))
        return float(value)
    if k in ("path", "str"):
        if not isinstance(value, str):
            raise ConfigError("%s takes a string, not %r" % (opt.name, value))
        return value
    if k == "enum":
        return _pick(opt, value, opt.choices, "value")
    if k == "floor" and isinstance(value, bool):
        return value
    if k in ("list", "floor", "bindings"):
        if isinstance(value, str) or not isinstance(value, (list, tuple)):
            raise ConfigError("%s takes a list, e.g. [\"A\", \"B\"], not %r"
                              % (opt.name, value)
                              + (" (or true / false)" if k == "floor"
                                 else ""))
        out = []
        for v in value:
            if not isinstance(v, str):
                raise ConfigError("%s: the items are strings, not %r"
                                  % (opt.name, v))
            allowed = _items(opt)
            if allowed and opt.open:
                v = next((a for a in allowed if a.upper() == v.upper()),
                         v.upper())
            elif allowed:
                v = _pick(opt, v, allowed, "item")
            if v not in out:
                out.append(v)
        return out
    raise AssertionError(k)


def defaults():
    """{section: {key: default}}"""
    out = {s: {} for s in SECTIONS}
    for o in OPTIONS:
        d = o.default
        out[o.section][o.key] = list(d) if isinstance(d, list) else d
    return out


def to_argv(values):
    """the generator's arguments for {section: {key: value}}: one flag per
    option that differs from its default (the internal bridge to the
    generator, docs/config-design.md section 4)"""
    argv = []
    for o in OPTIONS:
        v = values.get(o.section, {}).get(o.key, o.default)
        if v == o.default or v is None:
            continue
        if o.kind == "bindings":
            argv += ["--" + b for b in v]
        elif o.kind == "bool":
            argv.append(o.flag)
        elif o.kind == "notbool":
            argv.append(o.flag)          # true by default: only false here
        elif o.kind == "floor":
            argv.append(o.flag if v is True
                        else "%s=%s" % (o.flag, ",".join(v)))
        elif o.kind == "list":
            argv.append("%s=%s" % (o.flag, ",".join(v)))
        else:
            argv.append("%s=%s" % (o.flag, v))
    return argv


def from_namespace(ns):
    """{section: {key: value}} from the generator's parsed arguments (the
    inverse of to_argv; used to convert option files)"""
    out = defaults()
    for o in OPTIONS:
        if o.kind == "bindings":
            v = [b for b in ("fortran", "python") if getattr(ns, b, False)]
        else:
            v = getattr(ns, o.dest)
            if o.kind == "notbool":
                v = not v
            elif o.kind == "floor":
                v = (False if not v else True if v == "DEFAULT"
                     else [x.strip() for x in v.split(",") if x.strip()])
            elif o.kind == "list":
                v = [x.strip() for x in (v or "").split(",") if x.strip()]
            elif o.kind == "path":
                v = str(v) if v else None
        out[o.section][o.key] = v
    return out


# -- TOML ---------------------------------------------------------------------

SECTION_DOC = {
    "project": "the network, the project directory and what is generated",
    "rates": "how the rate coefficients and the rate equations are written",
    "thermal": "the thermal processes, fixed in the generated code",
    "radiation": "photon bins, shielding and column densities",
    "dust": "dust tables and the dust state",
    "surface": "grain-surface chemistry",
    "solver": "the project's default tolerances (the calling code can "
              "still change them at run time)",
}


def to_toml(values, meta=None, header=""):
    """the TOML text of {section: {key: value}}, every option written (a
    None value as a comment: it keeps its default), plus [meta] when
    given"""
    import tomlio
    secs = {}
    if meta:
        secs["meta"] = dict(meta)
    for s in SECTIONS:
        secs[s] = {o.key: values.get(s, {}).get(o.key, o.default)
                   for o in OPTIONS if o.section == s}
    docs = dict(SECTION_DOC, meta="written by KROME-X: the version and the "
                "network this project was generated from")
    text = tomlio.dumps(secs, docs)
    return (header + "\n\n" if header else "") + text


# descriptions of the options that differ from a generator flag (the
# positive form of a negative flag, the bindings list)
HELP = {
    "project.network": "the .knet network file, or a legacy react_* file "
    "(converted into the project); required (legacy -n)",
    "project.outdir": "the project directory (default generated/, or "
    "build_NAME for a shipped test)",
    "project.bindings": "the optional bindings over the C interface: "
    "fortran (module krome, example.f90) and python (krome.py, example.py, "
    "make shared)",
    "project.runtime": "copy the runtime headers and write the Makefile "
    "and the example into the project; false: the bare generated files "
    "only",
    "surface.competition": "two-body surface rates with reaction-diffusion "
    "competition; false: without it (legacy)",
    "surface.tunnelling": "two-body surface reactions with a barrier Ea by "
    "the faster of thermal hopping and tunnelling through a 1 A rectangular "
    "barrier; false: thermal hopping only, exp(-Ea/Tdust) (the legacy "
    "Semenov framework)",
}


def help_text(opt):
    """the description of an option (from the generator's argument help
    until the flags are removed, docs/config-design.md step 6)"""
    if opt.name in HELP:
        return HELP[opt.name]
    for a in _codegen().build_parser()._actions:
        if a.dest == opt.dest:
            return a.help or ""
    return ""


def _example(opt, examples):
    k, d = opt.kind, opt.default
    raw = examples.get(opt.dest)
    if k == "bindings":
        return ["fortran", "python"]
    if k in ("bool", "notbool"):
        return not d
    if k == "floor":
        return ["Z", "H2"]
    if k == "list":
        return [x for x in (raw or "A,B").split(",") if x]
    if k == "enum":
        return next(c for c in opt.choices if c != d)
    if k == "int":
        return int(raw) if raw else 10
    if k == "float":
        return float(raw) if raw else 1e-6
    return raw or "file.txt"


def template():
    """the commented template: every option, with its description, default
    and an example value, grouped by section"""
    import textwrap
    import tomlio
    cg = _codegen()
    examples = dict(cg.OPTION_EXAMPLES)
    out = ["# KROME-X configuration: every generation option, commented out.",
           "# Copy this file, uncomment the options you need (the values "
           "shown are",
           "# examples) and generate the project with",
           "#",
           "#     python3 -m kromex myrun.toml",
           "#",
           "# or load it in Python: kromex.Preprocessor.from_file("
           "\"myrun.toml\").",
           "# An option left commented keeps its default.",
           "#",
           "# Written by python3 -m kromex --template FILE (options_template.toml "
           "is the",
           "# shipped copy)."]
    for s in SECTIONS:
        out += ["", "# " + "=" * 72, "# " + SECTION_DOC[s], "# " + "=" * 72,
                "[%s]" % s]
        for o in OPTIONS:
            if o.section != s:
                continue
            out.append("")
            out += ["# " + ln for ln in textwrap.wrap(help_text(o), 74)]
            d = o.default
            out.append("# default: %s" % (
                "none" if d is None and o.key == "network" else
                "see the description" if d is None else tomlio.value(d)))
            out.append("#%s = %s" % (o.key,
                                     tomlio.value(_example(o, examples))))
    return "\n".join(out) + "\n"
