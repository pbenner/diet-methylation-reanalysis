#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ -x .tools/r/usr/lib/R/bin/R ]]; then
  exec .tools/r/usr/lib/R/bin/R --vanilla --slave --file="$1" --args "${@:2}"
fi
exec Rscript --vanilla "$@"
