"""The KROME-X command line (docs/config-design.md section 6):

    python3 -m kromex run.toml                       generate the project
    python3 -m kromex run.toml --set thermal.gamma=CONST --set project.outdir=build2
    python3 -m kromex --test collapse [--set ...]    a shipped test (examples/tests)
    python3 -m kromex --tests                        list the shipped tests
    python3 -m kromex --template FILE                the commented template
    python3 -m kromex --check run.toml               check the network only
    python3 -m kromex --version

--set SECTION.KEY=VALUE: VALUE in TOML syntax (1e-6, true, ["H2", "Z"]);
a bare word is a string (--set thermal.gamma=EXACT)."""

import argparse
import sys

if sys.version_info < (3, 8):
    sys.exit("KROME-X needs Python 3.8 or newer (this is %d.%d)"
             % sys.version_info[:2])

import codegen  # noqa: E402 (the package put its folder on the path)
import tomlio  # noqa: E402
from preprocessor import Preprocessor, ConfigError, GenerationError  # noqa


def parse_set(item):
    """SECTION.KEY=VALUE -> (name, value)"""
    if "=" not in item:
        raise ConfigError("--set %s: expected SECTION.KEY=VALUE" % item)
    name, raw = item.split("=", 1)
    try:
        value = tomlio.loads("x = " + raw)["x"]
    except tomlio.TOMLDecodeError:
        value = raw            # a bare word: a string
    return name.strip(), value


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="kromex", description="KROME-X: generate a project from a "
        "configuration file (TOML); docs/manual.md section 2",
        epilog="--set SECTION.KEY=VALUE overrides an option (VALUE in TOML "
               "syntax; a bare word is a string). The options: "
               "options_template.toml, or --template FILE.")
    ap.add_argument("config", nargs="?", help="the configuration (TOML)")
    ap.add_argument("--set", action="append", default=[],
                    metavar="SECTION.KEY=VALUE", help="override an option")
    ap.add_argument("--test", metavar="NAME",
                    help="set up the shipped test NAME (examples/tests)")
    ap.add_argument("--tests", action="store_true",
                    help="list the shipped tests")
    ap.add_argument("--template", metavar="FILE",
                    help="write the commented template of every option")
    ap.add_argument("--check", action="store_true",
                    help="check the network with these options, write "
                         "nothing (exit status 1 on a failure)")
    ap.add_argument("--version", action="version",
                    version="KROME-X %s" % codegen.krome_version())
    a = ap.parse_args(argv)
    try:
        if a.tests:
            print("\n".join(Preprocessor.tests()))
            return 0
        if a.template:
            with open(a.template, "w") as f:
                f.write(Preprocessor.template())
            print("wrote %s (the options template)" % a.template)
            return 0
        if a.test and a.config:
            ap.error("give a configuration file or --test, not both")
        if a.test:
            kp = Preprocessor.from_test(a.test)
        elif a.config:
            kp = Preprocessor.from_file(a.config)
        else:
            ap.error("no configuration file (or --test NAME)")
        for item in a.set:
            kp.set(*parse_set(item))
        if a.check:
            return 0 if kp.check() else 1
        kp.run()
        return 0
    except (ConfigError, GenerationError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
