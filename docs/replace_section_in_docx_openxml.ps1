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

Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

$wNs = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

function Get-ParagraphText {
    param(
        [System.Xml.XmlNode]$Paragraph,
        [System.Xml.XmlNamespaceManager]$Ns
    )
    $parts = New-Object System.Collections.Generic.List[string]
    foreach ($node in $Paragraph.SelectNodes('.//w:t | .//w:tab | .//w:br', $Ns)) {
        if ($node.LocalName -eq 't') {
            [void]$parts.Add($node.InnerText)
        } elseif ($node.LocalName -eq 'tab') {
            [void]$parts.Add("`t")
        } elseif ($node.LocalName -eq 'br') {
            [void]$parts.Add("`n")
        }
    }
    return (-join $parts)
}

function New-TextParagraph {
    param(
        [xml]$Doc,
        [string]$Text,
        [string]$StyleId,
        [int]$RevisionId,
        [string]$RevisionDate
    )

    $p = $Doc.CreateElement('w', 'p', $wNs)
    $pPr = $Doc.CreateElement('w', 'pPr', $wNs)
    $pStyle = $Doc.CreateElement('w', 'pStyle', $wNs)
    $pStyle.SetAttribute('val', $wNs, $StyleId)
    [void]$pPr.AppendChild($pStyle)
    [void]$p.AppendChild($pPr)

    $ins = $Doc.CreateElement('w', 'ins', $wNs)
    $ins.SetAttribute('id', $wNs, [string]$RevisionId)
    $ins.SetAttribute('author', $wNs, 'Codex')
    $ins.SetAttribute('date', $wNs, $RevisionDate)

    $r = $Doc.CreateElement('w', 'r', $wNs)
    $t = $Doc.CreateElement('w', 't', $wNs)
    $t.SetAttribute('space', 'http://www.w3.org/XML/1998/namespace', 'preserve')
    $t.InnerText = $Text
    [void]$r.AppendChild($t)
    [void]$ins.AppendChild($r)
    [void]$p.AppendChild($ins)
    return $p
}

function Add-ZipEntryFromBytes {
    param(
        [System.IO.Compression.ZipArchive]$Archive,
        [string]$Name,
        [byte[]]$Bytes
    )
    $entry = $Archive.CreateEntry($Name, [System.IO.Compression.CompressionLevel]::Optimal)
    $stream = $entry.Open()
    try {
        $stream.Write($Bytes, 0, $Bytes.Length)
    } finally {
        $stream.Dispose()
    }
}

function Read-ZipEntryText {
    param(
        [System.IO.Compression.ZipArchiveEntry]$Entry
    )
    $stream = $Entry.Open()
    try {
        $reader = New-Object System.IO.StreamReader($stream, [System.Text.Encoding]::UTF8)
        return $reader.ReadToEnd()
    } finally {
        $stream.Dispose()
    }
}

$startHeading = [string]::Concat(
    [char]0x7B2C, [char]0x4E00, [char]0x8282,
    ' Condensate-quant',
    [char]0x5206, [char]0x6790, [char]0x6D41, [char]0x7A0B
)
$endHeading = [string]::Concat(
    [char]0x7B2C, [char]0x4E8C, [char]0x8282,
    ' ',
    [char]0x8F6C, [char]0x5F55, [char]0x7206, [char]0x53D1,
    [char]0x52A8, [char]0x529B, [char]0x5B66,
    [char]0x5206, [char]0x6790
)

$sourceArchive = [System.IO.Compression.ZipFile]::OpenRead($SourceDocx)
try {
    $documentEntry = $sourceArchive.GetEntry('word/document.xml')
    if ($null -eq $documentEntry) {
        throw 'word/document.xml not found.'
    }

    [xml]$documentXml = Read-ZipEntryText $documentEntry
    $ns = New-Object System.Xml.XmlNamespaceManager($documentXml.NameTable)
    $ns.AddNamespace('w', $wNs)

    $paragraphs = @($documentXml.SelectNodes('//w:body/w:p', $ns))
    $startIndex = -1
    $endIndex = -1
    for ($i = 0; $i -lt $paragraphs.Count; $i++) {
        $text = (Get-ParagraphText $paragraphs[$i] $ns).Trim()
        if ($text -eq $startHeading) {
            $startIndex = $i
        }
    }
    if ($startIndex -lt 0) {
        throw 'Could not locate first section body heading.'
    }
    for ($i = $startIndex + 1; $i -lt $paragraphs.Count; $i++) {
        $text = (Get-ParagraphText $paragraphs[$i] $ns).Trim()
        if ($text -eq $endHeading) {
            $endIndex = $i
            break
        }
    }
    if ($endIndex -lt 0) {
        throw 'Could not locate second section body heading.'
    }

    $body = $documentXml.SelectSingleNode('//w:body', $ns)
    $beforeNode = $paragraphs[$endIndex]

    for ($i = $startIndex; $i -lt $endIndex; $i++) {
        [void]$body.RemoveChild($paragraphs[$i])
    }

    $lines = Get-Content -LiteralPath $UpdateText -Encoding UTF8
    $now = (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ')
    $revisionId = 1000
    foreach ($line in $lines) {
        $trimmed = $line.Trim()
        $styleId = 'Normal'
        if ($trimmed -match '^\d+\.\d+\.\d+\s') {
            $styleId = 'Heading3'
        } elseif ($trimmed -match '^\d+\.\d+\s') {
            $styleId = 'Heading2'
        } elseif ($revisionId -eq 1000) {
            $styleId = 'Heading1'
        }
        $newP = @(New-TextParagraph $documentXml $line $styleId $revisionId $now) | Where-Object { $_ -is [System.Xml.XmlNode] } | Select-Object -Last 1
        [void]$body.InsertBefore($newP, $beforeNode)
        $revisionId++
    }

    $settingsText = $null
    $settingsEntry = $sourceArchive.GetEntry('word/settings.xml')
    if ($null -ne $settingsEntry) {
        [xml]$settingsXml = Read-ZipEntryText $settingsEntry
        $settingsNs = New-Object System.Xml.XmlNamespaceManager($settingsXml.NameTable)
        $settingsNs.AddNamespace('w', $wNs)
        $settingsRoot = $settingsXml.SelectSingleNode('/w:settings', $settingsNs)
        if ($null -eq $settingsXml.SelectSingleNode('/w:settings/w:trackRevisions', $settingsNs)) {
            $track = $settingsXml.CreateElement('w', 'trackRevisions', $wNs)
            [void]$settingsRoot.AppendChild($track)
        }
        $settingsText = $settingsXml.OuterXml
    }

    if (Test-Path -LiteralPath $OutputDocx) {
        Remove-Item -LiteralPath $OutputDocx -Force
    }

    $utf8 = New-Object System.Text.UTF8Encoding($false)
    $newDocumentBytes = $utf8.GetBytes($documentXml.OuterXml)
    $newSettingsBytes = if ($null -ne $settingsText) { $utf8.GetBytes($settingsText) } else { $null }

    $outArchive = [System.IO.Compression.ZipFile]::Open($OutputDocx, [System.IO.Compression.ZipArchiveMode]::Create)
    try {
        foreach ($entry in $sourceArchive.Entries) {
            if ($entry.FullName -eq 'word/document.xml') {
                Add-ZipEntryFromBytes $outArchive $entry.FullName $newDocumentBytes
            } elseif ($entry.FullName -eq 'word/settings.xml' -and $null -ne $newSettingsBytes) {
                Add-ZipEntryFromBytes $outArchive $entry.FullName $newSettingsBytes
            } else {
                $inStream = $entry.Open()
                try {
                    $newEntry = $outArchive.CreateEntry($entry.FullName, [System.IO.Compression.CompressionLevel]::Optimal)
                    $outStream = $newEntry.Open()
                    try {
                        $inStream.CopyTo($outStream)
                    } finally {
                        $outStream.Dispose()
                    }
                } finally {
                    $inStream.Dispose()
                }
            }
        }
    } finally {
        $outArchive.Dispose()
    }
} finally {
    $sourceArchive.Dispose()
}

Write-Output "Saved: $OutputDocx"



