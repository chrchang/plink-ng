#!/usr/bin/env python3
"""Rewrite a Bioconda plink/plink2 recipe for a plink-ng release tag.

BiocondaBot's autobump ignores any upstream version containing a "-", and
every plink-ng tag has one, so the bot never sees a new release.  This maps
the tag to the conda version the recipe uses and rewrites the recipe's
version, tag, sha256 and build number.

  bioconda_bump.py TAG META_YAML

Prints "recipe=..." and "version=..." lines (for $GITHUB_OUTPUT), plus
"bump=true" when META_YAML was rewritten, or "bump=false" when the release is
not newer than the recipe (e.g. an a.6 release after an a.7 one).  Must run
under a Python that has the conda package, for conda's own version ordering.
"""

import hashlib
import re
import sys
import urllib.request

from conda.models.version import VersionOrder

TARBALL_URL = "https://github.com/chrchang/plink-ng/archive/refs/tags/{}.tar.gz"


def tag_to_recipe_version(tag):
    """Maps a release tag to (recipe name, conda version).

    v2.0.0-a.7.8 -> plink2 2.0.0a.7.8   (the recipe's existing scheme)
    v1.9.0-rc3   -> plink  1.90rc3      (1.90b7.7 was v1.9.0-b.7.7, and
    v1.9.0-b.7.15 -> plink 1.90b7.15     "1.9.0rc3" would sort below it)
    """
    m = re.fullmatch(r"v2\.0\.0(?:-([a-z]+\.[0-9.]+))?", tag)
    if m:
        return "plink2", "2.0.0" + (m.group(1) or "")
    m = re.fullmatch(r"v1\.9\.([0-9]+)(?:-([a-z]+)\.?([0-9.]+))?", tag)
    if m:
        version = "1.9" + m.group(1)
        if m.group(2):
            version += m.group(2) + m.group(3)
        return "plink", version
    raise ValueError(f"unrecognized release tag {tag!r}")


def set_jinja_var(text, name, value):
    pattern = r'(\{%\s*set\s+' + name + r'\s*=\s*")[^"]*("\s*%\})'
    new_text, n = re.subn(pattern, lambda m: m.group(1) + value + m.group(2), text)
    if n != 1:
        raise ValueError(f"expected one '{{% set {name} = \"...\" %}}' line, found {n}")
    return new_text


def set_yaml_key(text, key, value):
    pattern = r"(^\s*" + key + r":\s*)\S+"
    new_text, n = re.subn(pattern, lambda m: m.group(1) + value, text, flags=re.M)
    if n != 1:
        raise ValueError(f"expected one '{key}:' line, found {n}")
    return new_text


def main():
    tag, meta_path = sys.argv[1:]
    recipe, version = tag_to_recipe_version(tag)
    print(f"recipe={recipe}")
    print(f"version={version}")

    with open(meta_path) as f:
        meta = f.read()
    m = re.search(r'\{%\s*set\s+version\s*=\s*"([^"]*)"', meta)
    if not m:
        raise ValueError(f"{meta_path} has no '{{% set version = \"...\" %}}' line")
    current = m.group(1)
    if VersionOrder(version) <= VersionOrder(current):
        print(f"{recipe} is already at {current}; {version} is not newer", file=sys.stderr)
        print("bump=false")
        return

    with urllib.request.urlopen(TARBALL_URL.format(tag)) as resp:
        sha256 = hashlib.sha256(resp.read()).hexdigest()

    meta = set_jinja_var(meta, "version", version)
    meta = set_jinja_var(meta, "tag", tag)
    meta = set_yaml_key(meta, "sha256", sha256)
    meta = set_yaml_key(meta, "number", "0")
    with open(meta_path, "w") as f:
        f.write(meta)
    print(f"{recipe}: {current} -> {version} (sha256 {sha256})", file=sys.stderr)
    print("bump=true")


if __name__ == "__main__":
    main()
