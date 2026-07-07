#!/usr/bin/env bash
# View local STEP files on a phone (or any tailnet device) via a local step-loupe
# build served over Tailscale.
#
# Serves a copy of ../step-loupe/step-loupe.html plus the given STEP files from a
# static HTTP server bound to this machine's Tailscale IP. The page and the file
# share one HTTP origin, so no CORS / HTTPS / mixed-content applies — just open
# the printed URL on a device connected to the same tailnet.
#
#   tools/loupe-serve.sh path/to/model.step [more.step ...]
#
# With no arguments it only ensures the server is running and prints its URL.
set -euo pipefail

PORT=8787
DIR="$HOME/.cache/nacre-loupe"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOUPE_SRC="$REPO_ROOT/../step-loupe/step-loupe.html"

mkdir -p "$DIR"

# Refresh the viewer, building the single-file bundle first if it is missing.
if [ ! -f "$LOUPE_SRC" ] && [ -f "$REPO_ROOT/../step-loupe/scripts/build_single.py" ]; then
    (cd "$REPO_ROOT/../step-loupe" && python3 scripts/build_single.py)
fi
if [ ! -f "$LOUPE_SRC" ]; then
    echo "error: step-loupe.html not found at $LOUPE_SRC (build it in ../step-loupe)" >&2
    exit 1
fi
cp "$LOUPE_SRC" "$DIR/step-loupe.html"

TSIP="$(tailscale ip -4 | head -1)"
if [ -z "$TSIP" ]; then
    echo "error: no Tailscale IPv4 address (is tailscale up?)" >&2
    exit 1
fi

# Start the server once; reuse it if the port already answers (idempotent).
if ! curl -sf -o /dev/null "http://$TSIP:$PORT/step-loupe.html" 2>/dev/null; then
    nohup python3 -m http.server "$PORT" --directory "$DIR" --bind "$TSIP" >/dev/null 2>&1 &
    disown 2>/dev/null || true
    curl -sf -o /dev/null --retry 10 --retry-connrefused --retry-delay 1 --retry-max-time 15 \
        "http://$TSIP:$PORT/step-loupe.html" \
        || { echo "error: server did not come up on $TSIP:$PORT" >&2; exit 1; }
fi

if [ "$#" -eq 0 ]; then
    echo "server up: http://$TSIP:$PORT/  (drop STEP files into $DIR or pass them as args)"
    exit 0
fi

for f in "$@"; do
    if [ ! -f "$f" ]; then
        echo "skip (not a file): $f" >&2
        continue
    fi
    name="$(basename "$f")"
    cp "$f" "$DIR/$name"
    echo "http://$TSIP:$PORT/step-loupe.html?file=$name"
done
