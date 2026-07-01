param(
    [string]$SourceDocx,
    [string]$UpdateText,
    [string]$OutputDocx
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $SourceDocx)) {
    throw "Source DOCX not found: $SourceDocx"
}
if (-not (Test-Path -LiteralPath $UpdateText)) {
    throw "Update text not found: $UpdateText"
}

Copy-Item -LiteralPath $SourceDocx -Destination $OutputDocx -Force

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$word.ScreenUpdating = $false

try {
    $doc = $word.Documents.Open($OutputDocx)
    $doc.TrackRevisions = $false

    $needleStart = [string]::Concat(
        [char]0x7B2C, [char]0x4E00, [char]0x8282,
        ' Condensate-quant',
        [char]0x5206, [char]0x6790, [char]0x6D41, [char]0x7A0B
    )
    $needleEnd = [string]::Concat(
        [char]0x7B2C, [char]0x4E8C, [char]0x8282,
        ' ',
        [char]0x8F6C, [char]0x5F55, [char]0x7206, [char]0x53D1,
        [char]0x52A8, [char]0x529B, [char]0x5B66,
        [char]0x5206, [char]0x6790
    )

    $search = $doc.Content.Duplicate()
    $starts = New-Object System.Collections.Generic.List[int]
    while ($search.Find.Execute($needleStart)) {
        [void]$starts.Add($search.Start)
        $search = $doc.Range($search.End, $doc.Content.End)
    }
    if ($starts.Count -lt 1) {
        throw "Could not find section start."
    }

    $startPos = $starts[$starts.Count - 1]
    $endSearch = $doc.Range($startPos, $doc.Content.End)
    if (-not $endSearch.Find.Execute($needleEnd)) {
        throw "Could not find section end."
    }
    $endPos = $endSearch.Start

    $replacement = Get-Content -LiteralPath $UpdateText -Raw -Encoding UTF8
    $replacement = $replacement -replace "`r`n", "`r"
    $replacement = $replacement -replace "`n", "`r"
    if (-not $replacement.EndsWith("`r")) {
        $replacement = $replacement + "`r"
    }

    $target = $doc.Range($startPos, $endPos)
    $target.Text = $replacement

    # Leave revision mode enabled for future manual edits.
    $doc.TrackRevisions = $true
    $doc.Save()
    $doc.Close()
} finally {
    try { $word.Quit() } catch {}
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}

Write-Output "Saved: $OutputDocx"
