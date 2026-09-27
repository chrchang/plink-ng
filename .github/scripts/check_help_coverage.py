#!/usr/bin/env python3
"""Check that every flag plink2 accepts is covered by plink2_help.cc.

Each HelpPrint() call in 2.0/plink2_help.cc starts with an "index": a
'\\0'-separated list of flag names for which `plink2 --help <name>` prints
that block.  This script extracts every flag name the command-line parser
recognizes and fails if any of them is missing from all index strings.

Only the index is checked.  A flag alias only needs to appear in the index;
it is fine for the printed help text to describe just one spelling.

Flag names are collected from:
  * 2.0/plink2_cmdline.cc: flags handled before the main parse (--help,
    --version, --script, --rerun, --d, --out, ...);
  * 2.0/plink2.cc, alias pass: `strequal_k(flagname_p, "<alias>", ...)`
    comparisons that rewrite an alias to its canonical name;
  * 2.0/plink2.cc, main parse: `strequal_k_unsafe(flagname_p2, "...")`
    comparisons inside `switch (*flagname_p)`, where flagname_p2 skips the
    first character (the case label), plus the StrStartsWithUnsafe()
    prefix families whose suffixes are matched with strcmp().

Retired flags, whose parser branch only prints an error (usually pointing to
the replacement), are not required to be in the index.

With --compare-website, the index is instead compared against the flag
search on cog-genomics.org, in both directions.

Standard library only; runs in well under a second, no build needed.
"""

import argparse
import os
import re
import sys

# Flags that are accepted by the parser but intentionally absent from the
# help index.  Every entry needs a reason.
ALLOWLIST = {
    # Conventional single-dash shortcuts handled in CmdlineParsePhase1: -h/-?
    # print the short usage message, -v/-V print the version.
    "h": "shortcut for the short usage message",
    "?": "shortcut for the short usage message",
    "v": "shortcut for --version",
    "V": "shortcut for --version",
    # --help itself: its syntax is the usage line printed with every help
    # request.
    "help": "the help system itself",
    # Accepted so PLINK 1.x scripts keep running, but it only prints a note
    # saying it no longer has any effect.
    "allow-no-sex": "PLINK 1.x no-op",
}

# If extraction finds fewer flags than this, the parser code has probably
# changed shape and the regexes below need updating; fail loudly instead of
# passing vacuously.
MIN_EXPECTED_FLAGS = 400

WEBSITE_BASE = "https://www.cog-genomics.org"

VALID_NAME_RE = re.compile(r"[0-9A-Za-z?][0-9A-Za-z_-]*")
HELPPRINT_RE = re.compile(r'HelpPrint\("((?:[^"\\]|\\.)*)"')
MAIN_SWITCH_RE = re.compile(r'^      switch \(\*flagname_p\) \{$')
CASE_RE = re.compile(r"^      case '(.)':")
DEFAULT_RE = re.compile(r'^      default:')
P2_EQUAL_RE = re.compile(r'strequal_k_unsafe\(flagname_p2, "([^"]*)"\)')
P2_PREFIX_RE = re.compile(r'StrStartsWithUnsafe\(flagname_p2, "([^"]*)"\)')
SUBFLAG_RE = re.compile(r'strcmp\(subflag, "([^"]*)"\)')
EMPTY_P2_RE = re.compile(r"\*flagname_p2 == '\\0'")
BRANCH_RE = re.compile(r"^        (?:\} else )?if \(")
BRANCH_END_RE = re.compile(r"^        \}")
RETIRED_ERROR_RE = re.compile(r'(?:logerrputs|logerrprintfww|logerrprintf|snprintf)\((?:g_logbuf, kLogbufSize, )?"Error: ')
RETIRED_GOTO_RE = re.compile(r"goto main_ret_INVALID_CMDLINE\w*;$")
ALIAS_RE = re.compile(r'strequal_k\(flagname_p, "([^"]*)", flag_slen\)')
CMDLINE_RES = (
    re.compile(r'strequal_k\(flagname_p, "([^"]*)", flagname_p_slen\)'),
    re.compile(r'strequal_k_unsafe\(cur_flag, "([^"]*)"\)'),
    re.compile(r'Memcmp\("([^"]*)", cur_flag, \d+\)'),
    re.compile(r'strcmp\("--([^"]*)", argvk\[arg_idx\]\)'),
)


def unescape_c(s):
    return s.encode("latin-1").decode("unicode_escape")


def help_index(help_src):
    names = set()
    for m in HELPPRINT_RE.finditer(help_src):
        for name in unescape_c(m.group(1)).split("\0"):
            if name:
                names.add(name)
    return names


def is_retired_body(body_lines):
    """True if a parser branch does nothing but print an error and bail out.

    That is how plink2 handles retired and not-yet-implemented flags (e.g.
    --freqx, --ibc, --merge-equal-pos): the flag name is recognized only so
    that the error can point somewhere useful.
    """
    saw_error = False
    for line in body_lines:
        s = line.strip()
        if not s or s.startswith("//"):
            continue
        if RETIRED_ERROR_RE.match(s):
            saw_error = True
        elif not RETIRED_GOTO_RE.match(s):
            return False
    return saw_error


def main_parse_flags(plink2_lines):
    """(accepted, retired) flag sets from the main switch in plink2.cc."""
    flags = set()
    retired = set()
    switch_starts = [i for i, line in enumerate(plink2_lines) if MAIN_SWITCH_RE.match(line)]
    if not switch_starts:
        sys.exit("error: could not locate the main flag switch in plink2.cc")
    # The alias pass also switches on *flagname_p (more deeply indented); the
    # main parse switch is the last one at this indentation.
    start = switch_starts[-1]
    letter = None
    prefix = None
    # Flags named in the condition of the current top-level if/else-if
    # branch of the current case, and that branch's body.
    branch_flags = set()
    branch_body = []
    in_condition = False

    def close_branch():
        if branch_flags and is_retired_body(branch_body):
            retired.update(branch_flags)

    for line in plink2_lines[start + 1:]:
        if DEFAULT_RE.match(line):
            close_branch()
            break
        m = CASE_RE.match(line)
        if m:
            close_branch()
            letter = m.group(1)
            prefix = None
            branch_flags, branch_body, in_condition = set(), [], False
            continue
        if letter is None:
            continue
        if BRANCH_RE.match(line) or BRANCH_END_RE.match(line):
            close_branch()
            branch_flags, branch_body = set(), []
            in_condition = bool(BRANCH_RE.match(line))
        elif not in_condition:
            branch_body.append(line)
        line_flags = set()
        for m in P2_EQUAL_RE.finditer(line):
            line_flags.add(letter + m.group(1))
            prefix = None
        m = P2_PREFIX_RE.search(line)
        if m:
            prefix = letter + m.group(1)
        if prefix is not None:
            for m in SUBFLAG_RE.finditer(line):
                line_flags.add(prefix + m.group(1))
        if EMPTY_P2_RE.search(line):
            line_flags.add(letter)
        flags |= line_flags
        if in_condition:
            branch_flags |= line_flags
            if line.rstrip().endswith("{"):
                in_condition = False
    return flags, retired


def website_index():
    """Flag names in the cog-genomics.org PLINK 2.0 flag search.

    The site's flag search box is driven by a static script (its filename
    carries a date stamp, so it is looked up from the front page) holding one
    entry per documented flag group, each with a `matches: ['flag', ...]`
    list.  That list plays the same role as a HelpPrint() index.
    """
    import urllib.request

    def fetch(url):
        req = urllib.request.Request(url, headers={"User-Agent": "plink-ng-ci"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read().decode("utf-8", "replace")

    front = fetch(WEBSITE_BASE + "/plink/2.0/")
    m = re.search(r'src="(/static/js/plink2_flag_search[^"]*\.js)"', front)
    if not m:
        sys.exit("error: flag search script not found on " + WEBSITE_BASE + "/plink/2.0/")
    script = fetch(WEBSITE_BASE + m.group(1))
    names = set()
    for m in re.finditer(r"matches: \[([^\]]*)\]", script):
        names.update(re.findall(r"'([^']*)'", m.group(1)))
    if len(names) < MIN_EXPECTED_FLAGS // 2:
        sys.exit(f"error: only {len(names)} flag names parsed from the website flag search; "
                 "its format has probably changed.")
    return names


def compare_website(help_src):
    index = help_index(help_src)
    web = website_index()
    only_web = sorted(web - index)
    only_help = sorted(index - web)
    print(f"Website flag search: {len(web)} names; plink2_help.cc index: {len(index)} names.")
    print("(Differences are expected when the posted binaries lag the repository head.)\n")
    print(f"{len(only_web)} name(s) on the website but in no HelpPrint() index:")
    for name in only_web:
        print(f"  --{name}")
    print(f"\n{len(only_help)} name(s) in a HelpPrint() index but not in the website flag search:")
    for name in only_help:
        print(f"  --{name}")
    return 1 if (only_web or only_help) else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--src", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "2.0"),
                        help="plink2 source directory (default: 2.0/ in this checkout)")
    parser.add_argument("--list", action="store_true",
                        help="print every extracted flag (retired ones marked) and exit")
    parser.add_argument("--compare-website", action="store_true",
                        help="instead, compare the HelpPrint() index against the flag search on "
                        "cog-genomics.org/plink/2.0/ (needs network access)")
    args = parser.parse_args()

    def read(name):
        with open(os.path.join(args.src, name), encoding="latin-1") as f:
            return f.read()

    plink2_src = read("plink2.cc")
    cmdline_src = read("plink2_cmdline.cc")
    help_src = read("plink2_help.cc")

    if args.compare_website:
        return compare_website(help_src)

    accepted, retired = main_parse_flags(plink2_src.splitlines())
    aliases = set(ALIAS_RE.findall(plink2_src))
    accepted |= aliases
    for regex in CMDLINE_RES:
        accepted |= set(regex.findall(cmdline_src))

    if args.list:
        for flag in sorted(accepted):
            print(flag + ("  (retired)" if flag in retired else ""))
        return 0

    if len(accepted) < MIN_EXPECTED_FLAGS:
        print(f"error: only {len(accepted)} flags extracted from the parser (expected at least "
              f"{MIN_EXPECTED_FLAGS}); the flag-parsing code has probably changed shape and "
              f"{os.path.basename(__file__)} needs updating.", file=sys.stderr)
        return 1

    index = help_index(help_src)
    required = accepted - retired
    missing = sorted(f for f in required if f not in index and f not in ALLOWLIST)
    stale_allow = sorted(f for f in ALLOWLIST if f not in required or f in index)
    # A typo such as "freqx\frqx\0" (\f instead of \0) silently turns two
    # index entries into one unmatchable name.
    malformed = sorted(n for n in index if not VALID_NAME_RE.fullmatch(n))

    status = 0
    if missing:
        print(f"{len(missing)} flag(s) accepted by plink2 but missing from every HelpPrint() index "
              "in 2.0/plink2_help.cc:")
        for flag in missing:
            kind = " (alias)" if flag in aliases else ""
            print(f"  --{flag}{kind}")
        print("\nAdd each name to the index string of the relevant HelpPrint() block (aliases only\n"
              "need to be in the index, not in the printed text), or write a new block.")
        status = 1
    if malformed:
        print("\nMalformed HelpPrint() index entries (check the \\0 separators):")
        for name in malformed:
            print(f"  {name!r}")
        status = 1
    if stale_allow:
        print("\nAllowlist entries in check_help_coverage.py that are no longer needed:")
        for flag in stale_allow:
            print(f"  --{flag}")
        status = 1
    if not status:
        print(f"OK: all {len(required)} accepted flags appear in the plink2_help.cc index "
              f"({len(ALLOWLIST)} allowlisted; {len(retired)} retired flags, which only print an "
              "error, not required).")
    return status


if __name__ == "__main__":
    sys.exit(main())
