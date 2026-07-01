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

try {
    $doc = $word.Documents.Open($OutputDocx)
    $doc.TrackRevisions = $true
    $doc.ShowRevisions = $true

    $needleStart = '第一节 Condensate-quant分析流程'
    $needleEnd = '第二节 转录爆发动力学分析'

    $search = $doc.Content.Duplicate()
    $starts = New-Object System.Collections.Generic.List[int]
    while ($search.Find.Execute($needleStart)) {
        [void]$starts.Add($search.Start)
        $nextStart = $search.End
        $search = $doc.Range($nextStart, $doc.Content.End)
    }

    if ($starts.Count -lt 1) {
        throw "Could not find section start: $needleStart"
    }

    # Use the last hit to skip the table-of-contents entry and target the body heading.
    $startPos = $starts[$starts.Count - 1]
    $endSearch = $doc.Range($startPos, $doc.Content.End)
    if (-not $endSearch.Find.Execute($needleEnd)) {
        throw "Could not find section end: $needleEnd"
    }
    $endPos = $endSearch.Start

    $target = $doc.Range($startPos, $endPos)
    $target.Select()
    $word.Selection.Delete()

    $heading1 = $doc.Styles.Item(-2)
    $heading2 = $doc.Styles.Item(-3)
    $heading3 = $doc.Styles.Item(-4)
    $normal = $doc.Styles.Item(-1)

    $lines = Get-Content -LiteralPath $UpdateText -Encoding UTF8
    foreach ($line in $lines) {
        $trimmed = $line.Trim()

        if ($trimmed -match '^第一节\s') {
            $word.Selection.Style = $heading1
        } elseif ($trimmed -match '^\d+\.\d+\s') {
            $word.Selection.Style = $heading2
        } elseif ($trimmed -match '^\d+\.\d+\.\d+\s') {
            $word.Selection.Style = $heading3
        } else {
            $word.Selection.Style = $normal
        }

        $word.Selection.TypeText($line)
        $word.Selection.TypeParagraph()
    }

    # Keep Track Changes enabled for the user to continue editing in revision mode.
    $doc.TrackRevisions = $true
    foreach ($toc in $doc.TablesOfContents) {
        try {
            $toc.Update()
        } catch {
            # Updating the table of contents is convenient but non-critical.
        }
    }

    $doc.Save()
    $doc.Close()
} finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}

Write-Output "Saved: $OutputDocx"
