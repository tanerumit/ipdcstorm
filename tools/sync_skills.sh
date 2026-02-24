#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# sync_skills.sh
# Sync a library of skills from:
#   /WS/_shared/skills/<skill>/skill.md
# into the current repo at:
#   ./tools/skills/<skill>/skill.md
#
# Usage:
#   ./scripts/sync_skills.sh                 # sync all skills found
#   ./scripts/sync_skills.sh r weathergenr   # sync only these skills
#
# Env overrides:
#   SKILLS_SRC_ROOT=/WS/_shared/skills
#   SKILLS_DST_ROOT=tools/skills
# ------------------------------------------------------------------------------

SKILLS_SRC_ROOT="${SKILLS_SRC_ROOT:-/WS/_shared/skills}"
SKILLS_DST_ROOT="${SKILLS_DST_ROOT:-tools/skills}"

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -d "$SKILLS_SRC_ROOT" ]] || die "Source root not found: $SKILLS_SRC_ROOT"

mkdir -p "$SKILLS_DST_ROOT"

# Build skill list:
# - If args provided: treat as skill directory names
# - Else: discover all directories under source root that contain skill.md
skills=()
if [[ $# -gt 0 ]]; then
  skills=("$@")
else
  while IFS= read -r -d '' f; do
    # f = /WS/_shared/skills/<skill>/skill.md
    skill_dir="$(dirname "$f")"
    skill_name="$(basename "$skill_dir")"
    skills+=("$skill_name")
  done < <(find "$SKILLS_SRC_ROOT" -mindepth 2 -maxdepth 2 -type f -name 'skill.md' -print0)
fi

[[ ${#skills[@]} -gt 0 ]] || die "No skills found (or none specified)."

synced=0
skipped=0

for skill in "${skills[@]}"; do
  src="$SKILLS_SRC_ROOT/$skill/skill.md"
  dst="$SKILLS_DST_ROOT/$skill/skill.md"

  if [[ ! -f "$src" ]]; then
    echo "SKIP: missing source $src"
    skipped=$((skipped + 1))
    continue
  fi

  mkdir -p "$(dirname "$dst")"
  cp "$src" "$dst"
  echo "OK  : $src -> $dst"
  synced=$((synced + 1))
done

echo "Done. Synced: $synced, Skipped: $skipped"