param(
    [string]$InputText,
    [string]$OutputDocx
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $InputText)) {
    throw "Input text not found: $InputText"
}

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$word.ScreenUpdating = $false

try {
    $doc = $word.Documents.Add()
    $text = Get-Content -LiteralPath $InputText -Raw -Encoding UTF8
    $text = $text -replace "`r`n", "`r"
    $text = $text -replace "`n", "`r"

    $range = $doc.Range()
    $range.Text = $text

    $doc.TrackRevisions = $true
    $doc.SaveAs2($OutputDocx, 16)
    $doc.Close()
} finally {
    try { $word.Quit() } catch {}
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}

Write-Output "Saved: $OutputDocx"
