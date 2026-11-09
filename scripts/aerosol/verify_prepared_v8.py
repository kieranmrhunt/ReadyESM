"""v8 (9 Nov 2026): for prepare_uoexeter_v7.py (window December 1992 to January 1995). Only the script pin, the
expected verdict strings and the labels change. The checks are those of v7. The v7 text follows.

v7 (Paul's note on 65718246): every path argument is resolved to an absolute path before use, so the fresh
re-preparation (run in WORK_DIR) cannot misread relative inputs; otherwise identical to v6.

v6: for prepare_uoexeter_v6.py (runtime fingerprint, Paul): additionally requires the running interpreter's
runtime_fingerprint() (computed by the pinned v6 script's own function, in this process) to equal the pinned digest,
the receipt's runtime record and the prepared file's runtime_fingerprint attribute exactly.

v5: for prepare_uoexeter_v5.py (deterministic numerics, Paul's 65711449 finding); the exact fresh-output
comparison is unchanged and required. Labels name the version.

v4: for prepare_uoexeter_v4.py (cross-file bounds identity, Paul's review of v3); otherwise as v3.

Bind the raw UOEXETER files -> receipt -> prepared file (Claude Code, 6 Nov 2026; Paul's review of PA2-C v2).

Fail closed unless:
  1. the three raw files (ext, ssa, asy), the solar file and the two RRTMGP band files hash to the pinned values below,
     and the receipt's `inputs` equal those actual hashes exactly (same names, same values);
  2. the receipt's script_sha256 equals the actual hash of the given prepare_uoexeter_v3.py, which equals the pin;
  3. the receipt's output_sha256 equals the actual hash of the prepared file, and its verdict is the v3 verdict;
  4. the prepared file's attributes input_sha256 (JSON) and preparation_script_sha256 equal the receipt's;
  5. a fresh run of prepare_uoexeter_v3.py from the raw files (a new process, python -I, into a new directory)
     reproduces every variable of the prepared file exactly (same names, dimensions, dtypes and values) and every
     global and variable attribute, and its own receipt equals the given one in every field except output_sha256
     (the HDF5 container bytes are reported, not required, to be identical).
usage: python3 -I verify_prepared_v8.py RAW_DIR SOLAR RRTMGP_DIR PREPARED.nc RECEIPT.json PREPARE_V6.py WORK_DIR
"""
import hashlib
import json
import os
import runpy
import subprocess
import sys

import netCDF4
import numpy as np

RAW = {
    "ext_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc":
        "2a3980ff420487f94816f53f3fa2b6518bd28104d30ad66d55dd5f24baabc422",
    "ssa_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc":
        "c7406b8651a629eaec6cf9d84614e6b838ba5ff0d0b948dc60bdc7000b6f5abe",
    "asy_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-2-2-1_gnz_175001-202312.nc":
        "1f3f3f984b5c788cabe8df682f7ff31205688c4733575ee198339976e9154f60",
}
SOLAR = "4b844a3324b74be16faf7272b64ed7853ea2c606c4db841a318366ad223ebcac"
RRTMGP = {"rrtmgp-gas-sw-g224.nc": "584f1dd41ea9fc07d4ee3754eb1dafbd46ad3161cd6fd20fa06b6922b6f0702e",
          "rrtmgp-gas-lw-g256.nc": "4048360199d1917ed8f2ccaae2ec097d0f990da3bbad9830337b739b4fa01be7"}
PREPARE_V3 = "635050f03360583544d9c26a8201866f20be6dc85f3c241205555819674ae499"
RUNTIME_DIGEST = "a8c6f758970951cc82bbbab1b910ede9e44a5fa809e1aca05782fedf074f1a03"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def contents(path):
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        vars_ = {k: (v.dimensions, str(v.dtype), {a: v.getncattr(a) for a in v.ncattrs()}, v[:].copy())
                 for k, v in ds.variables.items()}
    return attrs, dims, vars_


def same_attr(a, b):
    return np.array_equal(a, b) if isinstance(a, np.ndarray) or isinstance(b, np.ndarray) else a == b


def main(raw_dir, solar, rrtmgp, prepared, receipt_path, script, work):
    raw_dir, solar, rrtmgp, prepared, receipt_path, script, work = (
        os.path.abspath(p) for p in (raw_dir, solar, rrtmgp, prepared, receipt_path, script, work))
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        print(f"{'ok  ' if cond else 'FAIL'} {name} {detail}")
        ok = ok and bool(cond)

    actual = {n: sha256(os.path.join(raw_dir, n)) for n in RAW}
    check("raw ext, ssa, asy hashes", actual == RAW, str({n[:3]: h[:12] for n, h in actual.items()}))
    actual[os.path.basename(solar)] = sha256(solar)
    check("solar hash", actual[os.path.basename(solar)] == SOLAR)
    for n in RRTMGP:
        actual[n] = sha256(os.path.join(rrtmgp, n))
    check("rrtmgp hashes", all(actual[n] == h for n, h in RRTMGP.items()))
    receipt = json.load(open(receipt_path))
    check("receipt inputs == actual input hashes", receipt["inputs"] == actual,
          str(sorted(set(receipt["inputs"]) ^ set(actual))))
    script_sha = sha256(script)
    check("receipt script == actual v6 script == pin", receipt["script_sha256"] == script_sha == PREPARE_V3, script_sha[:12])
    out_sha = sha256(prepared)
    check("receipt output == actual prepared file", receipt["output_sha256"] == out_sha, out_sha[:12])
    check("receipt verdict", receipt.get("verdict") == "UOEXETER_PREPARE_V7_DONE")
    current = runpy.run_path(script, run_name="verify")["runtime_fingerprint"]() if script_sha == PREPARE_V3 else {}
    check("running runtime fingerprint == pinned digest", current.get("digest") == RUNTIME_DIGEST,
          str(current.get("digest", ""))[:12])
    check("receipt runtime == running runtime", receipt.get("runtime") == current)
    attrs, dims, vars_ = contents(prepared)
    check("prepared input_sha256 attribute == receipt", json.loads(attrs["input_sha256"]) == receipt["inputs"])
    check("prepared preparation_script_sha256 attribute == receipt",
          attrs["preparation_script_sha256"] == receipt["script_sha256"])
    check("prepared runtime_fingerprint attribute == receipt runtime",
          json.loads(attrs.get("runtime_fingerprint", "{}")) == receipt.get("runtime"))
    os.makedirs(work, exist_ok=False)
    again, again_receipt = os.path.join(work, "reprepared.nc"), os.path.join(work, "reprepared_receipt.json")
    run = subprocess.run([sys.executable, "-I", script, raw_dir, solar, rrtmgp, again, again_receipt],
                         cwd=work, capture_output=True, text=True)
    check("fresh preparation ran", run.returncode == 0 and "UOEXETER_PREPARE_V7_DONE" in run.stdout,
          (run.stdout + run.stderr)[-300:])
    if run.returncode == 0:
        a2, d2, v2 = contents(again)
        check("dimensions identical", dims == d2)
        check("global attributes identical", attrs.keys() == a2.keys() and all(same_attr(attrs[k], a2[k]) for k in attrs))
        check("variable set identical", vars_.keys() == v2.keys())
        diff = [k for k in vars_ if k in v2 and not (vars_[k][0] == v2[k][0] and vars_[k][1] == v2[k][1] and
                vars_[k][2].keys() == v2[k][2].keys() and all(same_attr(vars_[k][2][a], v2[k][2][a]) for a in vars_[k][2])
                and np.array_equal(vars_[k][3], v2[k][3]))]
        check("every variable identical (dims, dtype, attributes, values)", not diff, str(diff))
        r2 = json.load(open(again_receipt))
        rd = sorted(k for k in set(receipt) | set(r2) if k != "output_sha256" and receipt.get(k) != r2.get(k))
        check("fresh receipt == given receipt except output_sha256", not rd, str(rd))
        print(f"note container_bytes_identical={sha256(again) == out_sha}")
    print("VERIFY_PREPARED_V8_PASS" if ok else "VERIFY_PREPARED_V8_FAIL")
    return ok


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    sys.exit(0 if main(*sys.argv[1:]) else 1)
