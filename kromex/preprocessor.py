"""kromex.Preprocessor: the configuration of a KROME-X project, and its
generation (docs/config-design.md).

    from kromex import Preprocessor

    kp = Preprocessor(network="networks/primordialZ.knet")
    kp.thermal.cooling = ["H2", "Z", "CHEM"]
    kp.thermal.metals = ["C", "C+", "O"]
    kp.solver.atol = 1e-20
    kp.run("build_Z")
    net = kp.build()        # make shared, and the loaded Network

The options are grouped in sections (kp.project, kp.rates, kp.thermal,
kp.radiation, kp.dust, kp.surface, kp.solver); every value is checked when
it is set, and an unknown name is an error."""

import importlib.util
import subprocess
from pathlib import Path

import codegen
import config
import tomlio
from config import ConfigError

GenerationError = codegen.GenerationError
import paths  # noqa: E402

ROOT = paths.ROOT


class Section:
    """the options of one section, as attributes"""

    def __init__(self, name, values):
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_values", values)

    def __getattr__(self, key):
        if key.startswith("_"):
            raise AttributeError(key)
        config.lookup("%s.%s" % (self._name, key))
        return self._values[key]

    def __setattr__(self, key, value):
        opt = config.lookup("%s.%s" % (self._name, key))
        self._values[key] = config.check(opt, value)

    def __dir__(self):
        return config.section_keys(self._name)

    def __repr__(self):
        return "[%s]\n%s" % (self._name, "\n".join(
            "%s = %r" % (k, self._values[k])
            for k in config.section_keys(self._name)))


class Preprocessor:
    """a KROME-X project configuration: the network and the options of
    every section, checked when set"""

    def __init__(self, network=None):
        values = config.defaults()
        object.__setattr__(self, "_values", values)
        object.__setattr__(self, "_sections", {
            s: Section(s, values[s]) for s in config.SECTIONS})
        object.__setattr__(self, "_test", None)
        object.__setattr__(self, "_outdir", None)
        if network is not None:
            self.project.network = str(network)

    # -- the options -------------------------------------------------------

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._sections:
            return self._sections[name]
        raise ConfigError(_where(name))

    def __setattr__(self, name, value):
        if name in config.SECTIONS:
            raise ConfigError("set the options of a section one by one, "
                              "e.g. kp.%s.%s = ..." % (
                                  name, config.section_keys(name)[0]))
        raise ConfigError(_where(name))

    def __dir__(self):
        return list(config.SECTIONS) + [
            "set", "get", "help", "to_dict", "to_toml", "save", "check",
            "run", "build", "from_dict", "from_file", "from_test", "tests",
            "template"]

    def set(self, name, value):
        """kp.set("thermal.gamma", "EXACT"): the same as the attribute"""
        opt = config.lookup(name)
        self._values[opt.section][opt.key] = config.check(opt, value)

    def get(self, name):
        opt = config.lookup(name)
        return self._values[opt.section][opt.key]

    @staticmethod
    def help(name):
        """the description of the option SECTION.KEY"""
        opt = config.lookup(name)
        print("%s (default %r): %s" % (opt.name, opt.default,
                                       config.help_text(opt)))

    def to_dict(self):
        """{section: {key: value}}, a copy"""
        return {s: dict(v) for s, v in self._values.items()}

    @classmethod
    def from_dict(cls, d):
        """a Preprocessor from {section: {key: value}}"""
        kp = cls()
        for s, keys in d.items():
            if not isinstance(keys, dict):
                raise ConfigError("%r is not a section of options" % s)
            for k, v in keys.items():
                kp.set("%s.%s" % (s, k), v)
        return kp

    @classmethod
    def from_file(cls, path, base=None):
        """a Preprocessor from a TOML file. Relative paths in it (the
        network, coolant and tolerance files) are relative to base, by
        default the file's folder. A project record (krome.toml, with
        [meta]) regenerates its own project: its folder is the default
        outdir, and its test is set up again"""
        path = Path(path)
        try:
            d = tomlio.load(path)
        except tomlio.TOMLDecodeError as exc:
            raise ConfigError("%s: %s" % (path, exc)) from exc
        except OSError as exc:
            raise ConfigError("cannot read %s: %s" % (path, exc)) from exc
        meta = d.pop("meta", None)
        try:
            kp = cls.from_dict(d)
        except ConfigError as exc:
            raise ConfigError("%s: %s" % (path, exc)) from exc
        base = Path(base) if base is not None else path.resolve().parent
        v = kp._values
        for sec, key in (("project", "network"), ("solver", "atol_file"),
                         ("solver", "rtol_file")):
            if v[sec][key] and not Path(v[sec][key]).is_absolute():
                v[sec][key] = str(base / v[sec][key])
        v["thermal"]["coolant_files"] = [
            f if Path(f).is_absolute() else str(base / f)
            for f in v["thermal"]["coolant_files"]]
        if meta is not None:
            if v["project"]["outdir"] is None:
                v["project"]["outdir"] = str(path.resolve().parent)
            if meta.get("test"):
                object.__setattr__(kp, "_test", meta["test"])
            net, want = v["project"]["network"], meta.get("network_sha256")
            if net and want and Path(net).exists():
                import hashlib
                got = hashlib.sha256(Path(net).read_bytes()).hexdigest()
                if got != want:
                    print("WARNING: %s has changed since this project was "
                          "generated (SHA-256 differs)" % net)
        return kp

    def to_toml(self):
        """this configuration as TOML text (every option)"""
        return config.to_toml(self._values)

    def save(self, path):
        """write this configuration to a TOML file"""
        Path(path).write_text(self.to_toml())

    @staticmethod
    def template():
        """the commented template of every option (options_template.toml)"""
        return config.template()

    def __repr__(self):
        d0 = config.defaults()
        lines = ["%s.%s = %r" % (s, k, v)
                 for s in config.SECTIONS for k, v in self._values[s].items()
                 if v != d0[s][k]]
        head = "Preprocessor(%s)" % (
            "test %r" % self._test if self._test else "")
        return "\n  ".join([head] + lines)

    # -- the shipped tests -------------------------------------------------

    @staticmethod
    def tests():
        """the names of the shipped tests (examples/tests)"""
        return codegen.list_tests()

    @classmethod
    def from_test(cls, name):
        """the options and network of the shipped test NAME
        (examples/tests/NAME/options.toml, its paths relative to the KROME-X
        folder); run() then also writes its driver, plot script and data
        (into build_NAME by default)"""
        if name not in codegen.list_tests():
            raise ConfigError("unknown test %r (tests: %s)"
                              % (name, ", ".join(codegen.list_tests())))
        kp = cls.from_file(codegen.TESTS_DIR / name / "options.toml",
                           base=ROOT)
        if codegen.test_driver(name) == "f90" and \
                "fortran" not in kp.project.bindings:
            kp.project.bindings = kp.project.bindings + ["fortran"]
        object.__setattr__(kp, "_test", name)
        return kp

    # -- generation --------------------------------------------------------

    def _args(self, outdir=None, check=False):
        ap = codegen.build_parser()
        v = self.to_dict()
        net = v["project"]["network"]
        if net:  # a shipped network may be named: "primordial"
            v["project"]["network"] = str(paths.find(net, "networks",
                                                     suffix=".knet"))
        for k in ("atol_file", "rtol_file"):
            if v["solver"][k]:
                v["solver"][k] = str(paths.find(v["solver"][k]))
        argv = config.to_argv(v)
        if self._test:
            argv.append("--test=" + self._test)
        if outdir is not None:
            argv.append("--outdir=" + str(outdir))
        if check:
            argv.append("--check")
        return ap, ap.parse_args(argv)

    def check(self):
        """check the network with these options, writing nothing; returns
        True when every check passes (the report is printed)"""
        ap, args = self._args(check=True)
        return bool(codegen.generate(args, ap)["check"]["ok"])

    def run(self, outdir=None):
        """generate the project (into outdir, project.outdir, build_NAME
        for a test, or generated/); returns its directory"""
        ap, args = self._args(outdir or self._values["project"]["outdir"])
        info = codegen.generate(args, ap)
        object.__setattr__(self, "_outdir", Path(info["outdir"]))
        return self._outdir

    def build(self, make_args=()):
        """make the shared library of the generated project (run() first;
        the Python binding is added if the project has none) and return
        its loaded Network"""
        if self._outdir is None:
            raise GenerationError("no project yet: call run() first")
        if not (self._outdir / "krome.py").exists():
            b = self._values["project"]["bindings"]
            self._values["project"]["bindings"] = b + ["python"]
            self.run(self._outdir)
        r = subprocess.run(["make", "shared", *make_args], cwd=self._outdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise GenerationError("make shared failed in %s:\n%s"
                                  % (self._outdir,
                                     (r.stdout + r.stderr)[-3000:]))
        lib = next(self._outdir.glob("libkrome.*"))
        spec = importlib.util.spec_from_file_location(
            "krome_%s" % self._outdir.name, self._outdir / "krome.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.Network(lib)


def _where(name):
    """the error for an option set or read without its section"""
    hits = [o.name for o in config.OPTIONS if o.key == name]
    if hits:
        return ("%r is in a section: use %s" % (
            name, " or ".join("kp." + h for h in hits)))
    return "%r is not a section (sections: %s)" % (
        name, ", ".join(config.SECTIONS))

