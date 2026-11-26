"""Verify the integrity of archived experimental evidence."""

from pathlib import Path
import hashlib


MANIFESTS = (
    "protocols/target_confirmation_SHA256SUMS.txt",
    "protocols/ecological_bootstrap_SHA256SUMS.txt",
    "protocols/geometry_validation_SHA256SUMS.txt",
)


def verify_archives(root):
    root = Path(root).resolve()
    verified = {}

    for manifest in MANIFESTS:
        manifest_path = root / manifest

        for line in manifest_path.read_text().splitlines():
            if not line.strip():
                continue

            expected, relative = line.split(maxsplit=1)
            relative = relative.lstrip("*")
            path = (root / relative).resolve()

            if not path.is_relative_to(root):
                raise ValueError(
                    f"Path outside repository: {relative}"
                )

            actual = hashlib.sha256(
                path.read_bytes()
            ).hexdigest()

            if actual != expected:
                raise ValueError(
                    f"Checksum mismatch: {relative}"
                )

            verified[relative] = actual

    return verified


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    verified = verify_archives(root)

    for name in verified:
        print(f"{name}: OK")

    print(f"Verified {len(verified)} archived files.")
