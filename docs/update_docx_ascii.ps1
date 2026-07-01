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
    $sectionOnePrefix = [string]::Concat([char]0x7B2C, [char]0x4E00, [char]0x8282)

    $search = $doc.Content.Duplicate()
    $starts = New-Object System.Collections.Generic.List[int]
    while ($search.Find.Execute($needleStart)) {
        [void]$starts.Add($search.Start)
        $nextStart = $search.End
        $search = $doc.Range($nextStart, $doc.Content.End)
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

        if ($trimmed.StartsWith($sectionOnePrefix)) {
            $word.Selection.Style = $heading1
        } elseif ($trimmed -match '^\d+\.\d+\.\d+\s') {
            $word.Selection.Style = $heading3
        } elseif ($trimmed -match '^\d+\.\d+\s') {
            $word.Selection.Style = $heading2
        } else {
            $word.Selection.Style = $normal
        }

        $word.Selection.TypeText($line)
        $word.Selection.TypeParagraph()
    }

    $doc.TrackRevisions = $true
    foreach ($toc in $doc.TablesOfContents) {
        try {
            $toc.Update()
        } catch {
        }
    }

    $doc.Save()
    $doc.Close()
} finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}

Write-Output "Saved: $OutputDocx"
