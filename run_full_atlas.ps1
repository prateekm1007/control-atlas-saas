# Control Atlas — Batch Execution Engine
# Regenerates all entries from 001 to N

$chimeraPath = "C:\Program Files\ChimeraX 1.11\bin\ChimeraX.exe"
$entriesRoot = "C:\Users\Administrator\control-atlas\entries"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "   CONTROL ATLAS — BATCH REGENERATION" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# Find all auto scripts
$scripts = Get-ChildItem -Path $entriesRoot -Recurse -Filter "*_auto.cxc"

foreach ($script in $scripts) {
    $entryName = $script.Directory.Name
    Write-Host "Rendering: $entryName" -NoNewline

    $cmd = "& `"$chimeraPath`" --script `"$($script.FullName)`""
    
    # Execute (Interactive fallback mode for safety)
    Invoke-Expression $cmd | Out-Null
    
    # Check if artifact exists (Simple heuristic)
    $pngs = Get-ChildItem -Path $script.Directory -Filter "*.png"
    if ($pngs.Count -gt 0) {
        Write-Host " [OK]" -ForegroundColor Green
    } else {
        Write-Host " [FAIL]" -ForegroundColor Red
    }
}

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "   ATLAS REGENERATION COMPLETE" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
