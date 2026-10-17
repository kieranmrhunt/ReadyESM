#!/usr/bin/env python3
"""Download the pinned January 1993 ECCO state from the public mirror."""

import hashlib
import os
import sys
import tempfile
import urllib.request
from pathlib import Path


BASE_URL = (
    "https://github.com/NumericalEarth/NumericalEarthArtifacts/"
    "releases/download/data-v1"
)
ASSETS = {
    "THETA_1993_01.nc": (
        103717993,
        "0b5e88bc1e865a8bdc566d858cda24333a967705eb731af2171e8a941fd318db",
    ),
    "SALT_1993_01.nc": (
        103717966,
        "8e8bdcfcefed4ef0e81ab4be3d6b1e2e9c6b99b3036dfcd0e4922ac78f596f99",
    ),
    "EVEL_1993_01.nc": (
        103718003,
        "d3d2927d3af6d5860ab821447bf2aef44e303e3d7e5a1988d098c0d6b5b7f68d",
    ),
    "NVEL_1993_01.nc": (
        103718005,
        "2ea1731bd8f5c653e0c89f9936d79865d47360df7a2ff2d044b077464301badf",
    ),
    "SSH_1993_01.nc": (
        2110218,
        "95b48394993f4dd8e87e383b35c70253a58c2ae1e81e3f1c0196c464b2c85488",
    ),
    "SIheff_1993_01.nc": (
        2110190,
        "36699cf11a99402782473c5fc051aa3e23dabadc46381cc7c03837c205353440",
    ),
    "SIarea_1993_01.nc": (
        2110183,
        "04518f541120851fc16484b065f5be520612ee17f6ec620cc3645bbd93308717",
    ),
}
RECEIPT_NAME = "ECCO4Monthly_1993_01_RECEIPT.txt"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def receipt_text():
    lines = [
        "schema=readiesm_ecco_monthly_bundle_v1",
        "dataset=ECCO_V4r4_interpolated_monthly_0p50deg",
        "source_month=1993-01",
        "coverage_start=1993-01-01T00:00:00",
        "coverage_end=1993-02-01T00:00:00",
        "source=https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/tag/data-v1",
        "layout=verified_payloads_downloaded_from_numericalearth_mirror",
        "variables=EVEL,NVEL,SALT,SIarea,SIheff,SSH,THETA",
    ]
    for filename in (
        "THETA_1993_01.nc",
        "SALT_1993_01.nc",
        "EVEL_1993_01.nc",
        "NVEL_1993_01.nc",
        "SSH_1993_01.nc",
        "SIheff_1993_01.nc",
        "SIarea_1993_01.nc",
    ):
        size, digest = ASSETS[filename]
        lines.append(f"{filename}_bytes={size}")
        lines.append(f"{filename}_sha256={digest}")
    return "\n".join(lines) + "\n"


def fetch(url, destination, expected_size, expected_sha256):
    destination = Path(destination)
    if destination.exists():
        if destination.is_file() and destination.stat().st_size == expected_size:
            if sha256(destination) == expected_sha256:
                return "reused"
        raise RuntimeError(f"existing ECCO asset differs: {destination}")
    temporary = destination.with_name(destination.name + ".tmp")
    if temporary.exists():
        raise RuntimeError(f"temporary ECCO asset already exists: {temporary}")
    request = urllib.request.Request(url, headers={"User-Agent": "ReadyESM/0.1"})
    with urllib.request.urlopen(request) as source, open(temporary, "xb") as target:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            target.write(block)
    if temporary.stat().st_size != expected_size:
        raise RuntimeError(
            f"ECCO asset has {temporary.stat().st_size} bytes, expected {expected_size}"
        )
    actual = sha256(temporary)
    if actual != expected_sha256:
        raise RuntimeError(
            f"ECCO asset SHA256 mismatch: expected {expected_sha256}, found {actual}"
        )
    os.rename(temporary, destination)
    return "downloaded"


def download(output_directory, base_url=BASE_URL, assets=ASSETS):
    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)
    allowed = set(assets) | {RECEIPT_NAME}
    unexpected = sorted(path.name for path in output.iterdir() if path.name not in allowed)
    if unexpected:
        raise RuntimeError(f"ECCO output directory has unexpected entries: {unexpected}")
    actions = {}
    for filename in sorted(assets):
        expected_size, expected_sha256 = assets[filename]
        actions[filename] = fetch(
            f"{base_url}/{filename}",
            output / filename,
            expected_size,
            expected_sha256,
        )

    receipt = receipt_text()
    receipt_path = output / RECEIPT_NAME
    if receipt_path.exists():
        if receipt_path.read_text(encoding="utf-8") != receipt:
            raise RuntimeError(f"existing ECCO receipt differs: {receipt_path}")
    else:
        receipt_path.write_text(receipt, encoding="utf-8")
    print(
        "ECCO_JANUARY_MIRROR_PASS "
        f"files={len(assets)} downloaded={sum(v == 'downloaded' for v in actions.values())} "
        f"reused={sum(v == 'reused' for v in actions.values())} output={output}"
    )
    return actions


def self_test():
    with tempfile.TemporaryDirectory() as root:
        source = Path(root) / "source"
        output = Path(root) / "output"
        source.mkdir()
        assets = {}
        for index, filename in enumerate(("a.nc", "b.nc")):
            payload = (f"asset {index}\n" * (index + 1)).encode()
            path = source / filename
            path.write_bytes(payload)
            assets[filename] = (len(payload), hashlib.sha256(payload).hexdigest())
        actions = download(output, source.as_uri(), assets)
        if set(actions.values()) != {"downloaded"}:
            raise AssertionError("synthetic assets were not downloaded")
        actions = download(output, source.as_uri(), assets)
        if set(actions.values()) != {"reused"}:
            raise AssertionError("synthetic assets were not reused")
        with open(output / "a.nc", "ab") as stream:
            stream.write(b"corruption")
        try:
            download(output, source.as_uri(), assets)
            raise AssertionError("corrupt cached asset was accepted")
        except RuntimeError as error:
            if "differs" not in str(error):
                raise
    print("ECCO_JANUARY_MIRROR_SELF_TEST_PASS")


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        self_test()
    elif len(sys.argv) == 2:
        download(sys.argv[1])
    else:
        raise SystemExit(
            "usage: download_ecco_january_mirror.py OUTPUT_DIRECTORY"
        )
