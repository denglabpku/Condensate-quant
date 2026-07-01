param(
    [string]$InputText,
    [string]$OutputDocx
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $InputText)) {
    throw "Input text not found: $InputText"
}

Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

function Escape-XmlText {
    param([string]$Text)
    if ($null -eq $Text) { return '' }
    return [System.Security.SecurityElement]::Escape($Text)
}

function Add-EntryText {
    param(
        [System.IO.Compression.ZipArchive]$Archive,
        [string]$Name,
        [string]$Text
    )
    $entry = $Archive.CreateEntry($Name, [System.IO.Compression.CompressionLevel]::Optimal)
    $stream = $entry.Open()
    try {
        $encoding = New-Object System.Text.UTF8Encoding($false)
        $bytes = $encoding.GetBytes($Text)
        $stream.Write($bytes, 0, $bytes.Length)
    } finally {
        $stream.Dispose()
    }
}

$contentTypes = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
</Types>
'@

$rels = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
'@

$docRels = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/settings" Target="settings.xml"/>
</Relationships>
'@

$styles = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal">
    <w:name w:val="Normal"/>
    <w:qFormat/>
    <w:pPr><w:spacing w:after="120" w:line="276" w:lineRule="auto"/></w:pPr>
    <w:rPr><w:rFonts w:ascii="Arial" w:eastAsia="Microsoft YaHei" w:hAnsi="Arial"/><w:sz w:val="21"/></w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading1">
    <w:name w:val="heading 1"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr><w:keepNext/><w:spacing w:before="480" w:after="180"/><w:outlineLvl w:val="0"/></w:pPr>
    <w:rPr><w:b/><w:rFonts w:ascii="Arial" w:eastAsia="Microsoft YaHei" w:hAnsi="Arial"/><w:sz w:val="32"/></w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading2">
    <w:name w:val="heading 2"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr><w:keepNext/><w:spacing w:before="300" w:after="120"/><w:outlineLvl w:val="1"/></w:pPr>
    <w:rPr><w:b/><w:rFonts w:ascii="Arial" w:eastAsia="Microsoft YaHei" w:hAnsi="Arial"/><w:sz w:val="26"/></w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading3">
    <w:name w:val="heading 3"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr><w:keepNext/><w:spacing w:before="220" w:after="100"/><w:outlineLvl w:val="2"/></w:pPr>
    <w:rPr><w:b/><w:rFonts w:ascii="Arial" w:eastAsia="Microsoft YaHei" w:hAnsi="Arial"/><w:sz w:val="23"/></w:rPr>
  </w:style>
</w:styles>
'@

$settings = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:settings xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:trackRevisions/>
</w:settings>
'@

$lines = Get-Content -LiteralPath $InputText -Encoding UTF8
$body = New-Object System.Text.StringBuilder
$lineNo = 0
foreach ($line in $lines) {
    $lineNo++
    $trimmed = $line.Trim()
    $styleId = 'Normal'
    if ($trimmed -match '^第[一二三四五六七八九十]+节\s') {
        $styleId = 'Heading1'
    } elseif ($trimmed -match '^\d+\.\d+\.\d+\s') {
        $styleId = 'Heading3'
    } elseif ($trimmed -match '^\d+\.\d+\s') {
        $styleId = 'Heading2'
    } elseif ($lineNo -eq 1) {
        $styleId = 'Heading1'
    }
    $escaped = Escape-XmlText $line
    [void]$body.AppendLine("    <w:p><w:pPr><w:pStyle w:val=`"$styleId`"/></w:pPr><w:r><w:t>$escaped</w:t></w:r></w:p>")
}

$document = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
$($body.ToString())
    <w:sectPr>
      <w:pgSz w:w="11906" w:h="16838"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"/>
    </w:sectPr>
  </w:body>
</w:document>
"@

if (Test-Path -LiteralPath $OutputDocx) {
    Remove-Item -LiteralPath $OutputDocx -Force
}

$archive = [System.IO.Compression.ZipFile]::Open($OutputDocx, [System.IO.Compression.ZipArchiveMode]::Create)
try {
    Add-EntryText $archive '[Content_Types].xml' $contentTypes
    Add-EntryText $archive '_rels/.rels' $rels
    Add-EntryText $archive 'word/_rels/document.xml.rels' $docRels
    Add-EntryText $archive 'word/document.xml' $document
    Add-EntryText $archive 'word/styles.xml' $styles
    Add-EntryText $archive 'word/settings.xml' $settings
} finally {
    $archive.Dispose()
}

Write-Output "Saved: $OutputDocx"
