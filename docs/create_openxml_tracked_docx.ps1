param(
    [string]$UpdateText,
    [string]$OutputDocx
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $UpdateText)) {
    throw "Update text not found: $UpdateText"
}

Add-Type -AssemblyName System.IO.Compression.FileSystem

function Escape-XmlText {
    param([string]$Text)
    if ($null -eq $Text) {
        return ''
    }
    return [System.Security.SecurityElement]::Escape($Text)
}

function Write-Utf8NoBom {
    param(
        [string]$Path,
        [string]$Text
    )
    $encoding = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($Path, $Text, $encoding)
}

$workRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("codex_docx_" + [System.Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $workRoot | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workRoot "_rels") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workRoot "word") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workRoot "word\_rels") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workRoot "docProps") | Out-Null

try {
    $contentTypes = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
'@

    $rels = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
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
  <w:revisionView w:markup="1" w:comments="1" w:insDel="1" w:formatting="1"/>
</w:settings>
'@

    $now = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
    $core = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>图像分析代码说明 更新版</dc:title>
  <dc:creator>Codex</dc:creator>
  <cp:lastModifiedBy>Codex</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">$now</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">$now</dcterms:modified>
</cp:coreProperties>
"@

    $app = @'
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Microsoft Word</Application>
</Properties>
'@

    $lines = Get-Content -LiteralPath $UpdateText -Encoding UTF8
    $body = New-Object System.Text.StringBuilder
    $revisionId = 1
    foreach ($line in $lines) {
        $trimmed = $line.Trim()
        $styleId = 'Normal'
        if ($trimmed -match '^\d+\.\d+\.\d+\s') {
            $styleId = 'Heading3'
        } elseif ($trimmed -match '^\d+\.\d+\s') {
            $styleId = 'Heading2'
        } elseif ($revisionId -eq 1) {
            $styleId = 'Heading1'
        }

        $escaped = Escape-XmlText $line
        [void]$body.AppendLine("    <w:p>")
        [void]$body.AppendLine("      <w:pPr><w:pStyle w:val=`"$styleId`"/></w:pPr>")
        [void]$body.AppendLine("      <w:ins w:id=`"$revisionId`" w:author=`"Codex`" w:date=`"$now`"><w:r><w:t xml:space=`"preserve`">$escaped</w:t></w:r></w:ins>")
        [void]$body.AppendLine("    </w:p>")
        $revisionId++
    }

    $document = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
$($body.ToString())
    <w:sectPr>
      <w:pgSz w:w="11906" w:h="16838"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>
    </w:sectPr>
  </w:body>
</w:document>
"@

    Write-Utf8NoBom (Join-Path $workRoot "[Content_Types].xml") $contentTypes
    Write-Utf8NoBom (Join-Path $workRoot "_rels\.rels") $rels
    Write-Utf8NoBom (Join-Path $workRoot "word\_rels\document.xml.rels") $docRels
    Write-Utf8NoBom (Join-Path $workRoot "word\styles.xml") $styles
    Write-Utf8NoBom (Join-Path $workRoot "word\settings.xml") $settings
    Write-Utf8NoBom (Join-Path $workRoot "word\document.xml") $document
    Write-Utf8NoBom (Join-Path $workRoot "docProps\core.xml") $core
    Write-Utf8NoBom (Join-Path $workRoot "docProps\app.xml") $app

    if (Test-Path -LiteralPath $OutputDocx) {
        Remove-Item -LiteralPath $OutputDocx -Force
    }
    [System.IO.Compression.ZipFile]::CreateFromDirectory($workRoot, $OutputDocx)
} finally {
    if (Test-Path -LiteralPath $workRoot) {
        Remove-Item -LiteralPath $workRoot -Recurse -Force
    }
}

Write-Output "Saved: $OutputDocx"
