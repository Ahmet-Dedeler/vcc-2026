#!/usr/bin/env bash
# Pull a built arm's package from S3, verify it against the checksum the worker wrote,
# and submit it. Kept separate from the build so a package can be fetched, inspected and
# held without spending a slot -- we get two per UTC day and rank tracks the LATEST entry,
# so an accidental submission actively demotes us.
set -Eeuo pipefail

ARM="${1:?usage: fetch_and_submit.sh <arm> [model-name]}"
NAME="${2:-}"
B="${BUCKET:-<your-s3-bucket>}"
DEST="${DEST:-submissions}"
VCC="${DEST}/${ARM}.vcc"

mkdir -p "$DEST"
if [[ ! -s "$VCC" ]]; then
  echo "--- fetching ${ARM} ---"
  aws s3 cp "s3://${B}/output/exp094/${ARM}/${ARM}.vcc" "$VCC"
fi
aws s3 cp "s3://${B}/output/exp094/${ARM}/${ARM}.vcc.sha256" "${VCC}.sha256.remote" --only-show-errors

want="$(awk '{print $1}' "${VCC}.sha256.remote")"
got="$(shasum -a 256 "$VCC" | awk '{print $1}')"
if [[ "$want" != "$got" ]]; then
  echo "CHECKSUM MISMATCH for ${ARM}: local ${got} != worker ${want}" >&2
  exit 65
fi
echo "checksum OK: ${got}"
# The build metadata is small and carries the audited config -- including experiment 099's
# measured `f` under effects.cd4.residual_frame. Optional: arms built before fastbuild.sh
# started uploading it have none.
aws s3 cp "s3://${B}/output/exp094/${ARM}/${ARM}.h5ad.meta.json" "${VCC}.meta.json" \
  --only-show-errors 2>/dev/null \
  && echo "metadata: ${VCC}.meta.json" \
  || echo "no build metadata at this arm's prefix"
ls -la "$VCC"

if [[ -z "$NAME" ]]; then
  echo
  echo "No model name given -- fetched and verified only, NOT submitted."
  echo "To submit:  bash $0 ${ARM} '<model name>'"
  exit 0
fi

echo "--- submitting ${ARM} as '${NAME}' ---"
vcc submit "$VCC" --model-name "$NAME" --wait 2>&1 | tee "${DEST}/${ARM}.submit.log"
