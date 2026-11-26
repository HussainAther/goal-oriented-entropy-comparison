"""Retrieve reference metadata from DOI registries."""

from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import quote
from datetime import datetime, timezone
import json
import re
import time


ROOT = Path(__file__).resolve().parents[1]

bib = (ROOT / "references.bib").read_text()

entries = re.split(r"(?=@\w+\{)", bib)

records = []

for entry in entries:
    key = re.search(r"@\w+\{([^,]+),", entry)
    doi = re.search(r"doi\s*=\s*\{([^}]+)\}", entry)

    if not key:
        continue

    record = {"key": key.group(1)}

    if not doi:
        record["status"] = "no_doi"
        records.append(record)
        continue

    identifier = doi.group(1)
    record["doi"] = identifier

    if identifier.startswith("10.5061/"):
        url = (
            "https://api.datacite.org/dois/"
            + quote(identifier, safe="")
        )
        registry = "DataCite"
    else:
        url = (
            "https://api.crossref.org/works/"
            + quote(identifier, safe="")
        )
        registry = "Crossref"

    record["registry"] = registry
    record["url"] = url

    try:
        request = Request(
            url,
            headers={
                "User-Agent": "EntropyComparator/1.0",
                "Accept": "application/json",
            },
        )

        with urlopen(request, timeout=30) as response:
            payload = json.load(response)

        record["status"] = "retrieved"
        record["metadata"] = (
            payload["data"]["attributes"]
            if registry == "DataCite"
            else payload["message"]
        )

    except Exception as exc:
        record["status"] = "failed"
        record["error"] = str(exc)

    print(
        record["key"],
        registry,
        record["status"],
        flush=True,
    )

    records.append(record)
    time.sleep(1)

result = {
    "retrieved_at_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "records": records,
}

output = ROOT / "literature/registry_metadata.json"

output.write_text(
    json.dumps(result, indent=2, ensure_ascii=False)
    + "\n"
)

print("\nTotal references:", len(records))
print(
    "Retrieved:",
    sum(r["status"] == "retrieved" for r in records),
)
print(
    "Failed:",
    sum(r["status"] == "failed" for r in records),
)
print("Output:", output)
