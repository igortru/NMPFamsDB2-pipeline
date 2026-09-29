#!/usr/bin/env bash
# Download the IMG/VR v4.1 full set (IMG_VR_2022-12-19_7) from the JGI Data Portal and check md5.
# usage: download_imgvr.sh OUTDIR [TOKEN_FILE]
#   TOKEN_FILE: file holding the JGI Data Portal session token, as copied from
#               avatar -> "Copy My Session Token" (default ~/.jgi_token, chmod 600)
# Uses the POST /download_files/ endpoint (documented in the JGI API tutorial), one file per
# request; the response is a zip, unpacked with Python zipfile (zip64-safe), then md5-checked.
# JAMO ids and md5 from the JGI files API (dataset Custom_MPI-IMG_VR). The v4.1 nucleotide file
# is byte-identical to the v4.0 copy (same md5), so the on-disk copy listed under
# IMG_VR_2022-09-20_6 is used.
set -euo pipefail

outdir=${1:?usage: download_imgvr.sh OUTDIR [TOKEN_FILE]}
token_file=${2:-$HOME/.jgi_token}
[ -r "$token_file" ] || { echo "token file not readable: $token_file" >&2; exit 1; }
mkdir -p "$outdir"

dataset=Custom_MPI-IMG_VR
# jamo_id  md5  file_name  size_bytes
files="
632a1b312de7c323533daabd e37e6a541e10a75a88577763f75fc9cd IMGVR_all_nucleotides.fna.gz 48319030883
63a22c0c3b5d0133c73fb095 fff611fc3983ce2ce8e9f9a4d671ae0e IMGVR_all_proteins.faa.gz 31477831962
63a22c0c3b5d0133c73fb097 658ce076cd26729892552ed5876ad17c IMGVR_all_Sequence_information.tsv 5047197637
63a22c0c3b5d0133c73fb09b ba749fd171bbf5821dcddeeff3db343b README.txt 4332
"

md5_of() { if command -v md5 >/dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

status=0
while read -r id sum name size; do
    [ -n "$id" ] || continue
    out="$outdir/$name"
    if [ -s "$out" ] && [ "$(md5_of "$out")" = "$sum" ]; then
        echo "ok (already present)  $name"; continue
    fi
    echo "downloading $name ($size bytes)"
    zip="$outdir/$name.download.zip"
    code=$(curl -sS -L --retry 5 --retry-delay 30 -X POST \
        -H "Authorization: $(cat "$token_file")" \
        -H "Content-Type: application/json" -H "accept: application/json" \
        -d "{\"ids\": {\"$dataset\": {\"file_ids\": [\"$id\"]}}, \"api_version\": \"2\"}" \
        -o "$zip" -w '%{http_code}' \
        "https://files-download.jgi.doe.gov/download_files/") || code=curl_error
    if [ "$code" != "200" ]; then
        echo "FAILED $name: HTTP $code $(head -c 300 "$zip" 2>/dev/null)" >&2
        [ "$code" = "409" ] && echo "  file is on tape: request a restore (request_archived_files), then re-run" >&2
        rm -f "$zip"; status=1; continue
    fi
    python3 - "$zip" "$name" "$out.part" <<'EOF' || { echo "FAILED $name: could not extract from zip (kept $zip)" >&2; status=1; continue; }
import shutil, sys, zipfile
z, name, dest = sys.argv[1:]
with zipfile.ZipFile(z) as zf:
    hits = [m for m in zf.namelist() if m.rsplit("/", 1)[-1] == name]
    if not hits:
        sys.exit(f"{name} not in zip; members: {zf.namelist()[:10]}")
    with zf.open(hits[0]) as src, open(dest, "wb") as dst:
        shutil.copyfileobj(src, dst, 16 << 20)
EOF
    got=$(md5_of "$out.part")
    if [ "$got" = "$sum" ]; then mv "$out.part" "$out"; rm -f "$zip"; echo "ok  $name  md5=$got"
    else echo "FAILED $name: md5 $got, expected $sum (kept $zip)" >&2; status=1; fi
done <<< "$files"
exit $status
