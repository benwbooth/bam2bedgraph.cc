#!/usr/bin/env bash
set -eu
(cd "$(dirname "$0")" && make -s)
exec "$(dirname "$0")/$(basename "${0%.*}")" "$@"
