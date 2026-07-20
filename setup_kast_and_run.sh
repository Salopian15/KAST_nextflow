#!/usr/bin/env bash
#
# Builds a Docker image containing the KAST 1.0.1 release binary, installs a
# shim so `kast` works transparently on macOS, then runs the HGT pipeline.
#
#   ./setup_kast_and_run.sh
#   KMER=5 ./setup_kast_and_run.sh          # override k-mer length
#   SKIP_RUN=1 ./setup_kast_and_run.sh      # build + install shim only
#
set -euo pipefail

RUNDIR="${RUNDIR:-$HOME/Documents/hgt_run}"
IMAGE="${IMAGE:-kast:1.0.1}"
SHIM="${SHIM:-$HOME/scripts/kast}"
KMER="${KMER:-3}"
SKIP_RUN="${SKIP_RUN:-0}"
KAST_URL="https://github.com/martinjvickers/KAST/releases/download/1.0.1/kast_1.0.1.tar.gz"

log()  { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
die()  { printf '\033[31mERROR:\033[0m %s\n' "$1" >&2; exit 1; }

# --- k-mer sanity: KAST allocates a dense 4^k array -------------------------
# RAM needed ~= 4^k * 4 bytes.  k=14 -> 1 GB, k=15 -> 4 GB, k=16 -> 16 GB,
# k=25 -> ~4 PB (KAST refuses to start).  KAST's own default is 3.
if [ "$KMER" -gt 15 ]; then
    die "KMER=$KMER needs ~$(( 4 ** (KMER - 15) * 4 )) GB+ of RAM. KAST will refuse to run.
       Use a value KAST can handle (its default is 3; <=14 needs <=1 GB):
         KMER=3 $0"
fi

log "Pre-flight"
command -v docker   >/dev/null 2>&1 || die "docker not found - install Docker Desktop for Mac"
docker info         >/dev/null 2>&1 || die "Docker daemon not running - start Docker Desktop"
command -v nextflow >/dev/null 2>&1 || die "nextflow not found in PATH"
[ -d "$RUNDIR" ]                    || die "run directory not found: $RUNDIR"
[ -f "$RUNDIR/main.nf" ]            || die "main.nf not found in $RUNDIR"
echo "  run dir : $RUNDIR"
echo "  image   : $IMAGE"
echo "  shim    : $SHIM"
echo "  kmer    : $KMER"

BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT

log "Downloading KAST 1.0.1 release binary"
curl -fSL --retry 3 -o "$BUILD/kast.tar.gz" "$KAST_URL"
tar xzf "$BUILD/kast.tar.gz" -C "$BUILD"
[ -f "$BUILD/kast" ] || die "expected 'kast' inside the release tarball"
chmod +x "$BUILD/kast"
file "$BUILD/kast" | sed 's/^/  /'

log "Building image (linux/amd64)"
cat > "$BUILD/Dockerfile" <<'DOCKERFILE'
FROM debian:bookworm-slim
# KAST 1.0.1 is a statically linked Linux x86-64 binary; no runtime deps needed.
COPY kast /usr/local/bin/kast
RUN chmod +x /usr/local/bin/kast
ENTRYPOINT []
DOCKERFILE
docker build --platform linux/amd64 -t "$IMAGE" "$BUILD"

log "Verifying image"
docker run --rm --platform linux/amd64 "$IMAGE" kast --version

log "Installing shim"
mkdir -p "$(dirname "$SHIM")"
# Back up an existing real binary (not a previously installed shell shim).
if [ -e "$SHIM" ] && ! head -c 2 "$SHIM" 2>/dev/null | grep -q '#!'; then
    mv "$SHIM" "$SHIM.linux-x86_64.bak"
    echo "  backed up original -> $SHIM.linux-x86_64.bak"
fi
cat > "$SHIM" <<EOF
#!/usr/bin/env bash
# Auto-generated: runs KAST 1.0.1 inside Docker (linux/amd64, emulated on arm64).
# \$RUNDIR is bind-mounted so Nextflow's staged symlinks resolve inside the container.
exec docker run --rm --platform linux/amd64 \\
  -v "$RUNDIR":"$RUNDIR" -w "\$PWD" \\
  "$IMAGE" kast "\$@"
EOF
chmod +x "$SHIM"
echo "  wrote $SHIM"

log "Shim test (from $RUNDIR)"
cd "$RUNDIR"
"$SHIM" --version

case ":$PATH:" in
  *":$(dirname "$SHIM"):"*) ;;
  *) echo "  NOTE: $(dirname "$SHIM") is not in PATH; add it or the pipeline won't find kast" ;;
esac

if [ "$SKIP_RUN" = "1" ]; then
    log "SKIP_RUN=1 - stopping before the pipeline"
    exit 0
fi

log "Running pipeline (kmer_length=$KMER)"
nextflow run main.nf -resume \
    --kmer_length "$KMER" \
    -with-report report.html \
    -with-trace
