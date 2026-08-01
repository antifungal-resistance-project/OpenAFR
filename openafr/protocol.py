"""Load and hash-verify the frozen run protocol (protocol.yaml).

protocol.yaml holds the numbers that define a run — the docking box + search
effort, and the gate's pass bar. They are the single easiest thing to silently
tune into a passing result, so this module pins the file's SHA-256 and the gate
re-checks it at run time. That is the same anti-cheat guarantee the
pre-registration gives, but machine-readable so the docking script and the gate
read one source of truth instead of hand-copied constants.

No PyYAML dependency: the file is a flat, self-authored, hash-frozen key/value
map, so a small strict reader (this module) is enough and keeps the pipeline
runnable with just rdkit + the docking toolkit — nothing beyond what
environment.yml already pins. If PyYAML is added there later, _parse() can be
swapped for yaml.safe_load().
"""
import hashlib
import pathlib

_ROOT = pathlib.Path(__file__).resolve().parent.parent
PROTOCOL = _ROOT / "protocol.yaml"

# SHA-256 of protocol.yaml at freeze time. Regenerate with:
#   python -m openafr.protocol --freeze
PROTOCOL_SHA = "6eca16b1be2e5e33ae42d4b1b7a63a0e9b891a5dd6cbfb19642dd3e9a149aba4"


def _scalar(tok):
    tok = tok.strip()
    if tok.startswith("[") and tok.endswith("]"):
        inner = tok[1:-1].strip()
        return [_scalar(x) for x in inner.split(",")] if inner else []
    for cast in (int, float):
        try:
            return cast(tok)
        except ValueError:
            pass
    return tok


def _parse(text):
    """Strict two-level reader: top-level 'section:' headers, indented 'key: value'."""
    data, section = {}, None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].rstrip()          # drop comments + trailing space
        if not line.strip():
            continue
        if line[0] not in " \t":                       # top-level section header
            key, sep, rest = line.partition(":")
            if not sep or rest.strip():
                raise ValueError(f"protocol.yaml: expected 'section:' header, got {line!r}")
            section = data[key.strip()] = {}
        else:                                          # indented key: value
            if section is None:
                raise ValueError(f"protocol.yaml: value {line!r} before any section")
            key, sep, rest = line.partition(":")
            if not sep:
                raise ValueError(f"protocol.yaml: expected 'key: value', got {line!r}")
            section[key.strip()] = _scalar(rest)
    return data


def _digest(path=PROTOCOL):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def load(path=PROTOCOL):
    """Parse protocol.yaml into {section: {key: value}}."""
    return _parse(pathlib.Path(path).read_text())


def verify(path=PROTOCOL):
    """Print whether the protocol is unmodified since freezing; return True if so.

    Soft check, matching the pre-registration: it never blocks a run, it makes a
    tampered protocol impossible to miss in the gate output.
    """
    got = _digest(path)
    if got == PROTOCOL_SHA:
        print(f"protocol verified unmodified (sha256 {got[:16]}...)")
        return True
    print("!! protocol.yaml HAS BEEN MODIFIED SINCE IT WAS FROZEN !!")
    print(f"   expected {PROTOCOL_SHA[:16]}...  got {got[:16]}...")
    print("   Any 'pass' under this protocol is not credible until re-frozen.")
    return False


# Map the docking section onto the variable names the docking shell script uses.
_SHELL_KEYS = [
    ("CX", "center_x"), ("CY", "center_y"), ("CZ", "center_z"),
    ("SIZE", "size"), ("EXHAUST", "exhaustiveness"),
    ("MODES", "num_modes"), ("SEED", "seed"),
]


def _emit_shell(section="docking"):
    vals = load()[section]
    print("; ".join(f"{var}={vals[key]}" for var, key in _SHELL_KEYS))


def main(argv=None):
    import json
    import sys
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else "--show"
    if cmd == "--shell":
        _emit_shell(argv[1] if len(argv) > 1 else "docking")
    elif cmd == "--freeze":
        print(_digest())
    elif cmd == "--verify":
        return 0 if verify() else 1
    else:
        print(json.dumps(load(), indent=2))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
