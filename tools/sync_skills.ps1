# tools/sync_skills.ps1
$ErrorActionPreference = "Stop"

# EDIT if needed
$skillsSrcRoot = "C:\Users\taner\WS\_shared\skills"   # or "C:\WS\_shared\skills"

$repoRoot = (Get-Location).Path
$skillsDstRoot = Join-Path $repoRoot ".agents\skills"

if (!(Test-Path -LiteralPath $skillsSrcRoot)) {
  throw "Source root not found: $skillsSrcRoot"
}

# Find all skills in shared library: /skills/<skill>/(skill.md|SKILL.md)
$skillFiles = Get-ChildItem -LiteralPath $skillsSrcRoot -Recurse -File |
  Where-Object { $_.Name -in @("skill.md", "SKILL.md") }

if ($skillFiles.Count -eq 0) {
  throw "No skill.md or SKILL.md files found under: $skillsSrcRoot"
}

New-Item -ItemType Directory -Force -Path $skillsDstRoot | Out-Null

$syncCount = 0
foreach ($f in $skillFiles) {
  # Expect: ...\skills\<skill>\skill.md
  $skillName = Split-Path (Split-Path $f.FullName -Parent) -Leaf

  $dstDir = Join-Path $skillsDstRoot $skillName
  $dst = Join-Path $dstDir "SKILL.md"

  New-Item -ItemType Directory -Force -Path $dstDir | Out-Null

  # BOM-safe copy: read raw, strip BOM, write UTF-8 (no BOM)
  $content = Get-Content -Raw -LiteralPath $f.FullName
  $content = $content -replace "^\uFEFF", ""
  $content = $content.TrimStart()  # ensure YAML starts at byte 0
  Set-Content -Encoding utf8 -LiteralPath $dst -Value $content

  Write-Host "OK  : $($f.FullName) -> $dst"
  $syncCount++
}

Write-Host "Done. Synced: $syncCount"