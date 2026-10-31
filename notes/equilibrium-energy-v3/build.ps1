# Build the full note, then each standalone chapter (pdflatex + bibtex).
# Chapters import equation and section numbers from the full document,
# so main is built first.
#
#   .\build.ps1
#   .\build.ps1 -OutDir build_tmp

param(
    [string]$OutDir = 'build'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$Root = $PSScriptRoot
$TexLive = 'C:\texlive\2024\bin\windows'
$Pdflatex = Join-Path $TexLive 'pdflatex.exe'
$Bibtex = Join-Path $TexLive 'bibtex.exe'

if (-not (Test-Path -LiteralPath $Pdflatex)) {
    throw "pdflatex not found: $Pdflatex"
}
if (-not (Test-Path -LiteralPath $Bibtex)) {
    throw "bibtex not found: $Bibtex"
}

$env:PATH = "$TexLive;$env:PATH"
if ([string]::IsNullOrEmpty($env:BIBINPUTS)) {
    $env:BIBINPUTS = '.;'
} elseif ($env:BIBINPUTS -notmatch '(^|;)\.(;|$)') {
    $env:BIBINPUTS = ".;$env:BIBINPUTS"
}

Set-Location -LiteralPath $Root
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$Docs = @(
    'main',
    'ch1_mag_base',
    'ch2_mag_beta',
    'ch3_mag_phase',
    'ch4_es_base',
    'ch5_es_beta',
    'ch6_compare',
    'appA_bm'
)

function Show-LogTail {
    param(
        [string]$Path,
        [int]$Lines = 40
    )
    if (Test-Path -LiteralPath $Path) {
        Write-Host "----- $(Split-Path -Leaf $Path) (last $Lines lines) -----"
        Get-Content -LiteralPath $Path -Tail $Lines | ForEach-Object { Write-Host $_ }
    }
}

function Invoke-Pdflatex {
    param([string]$Name)

    $log = Join-Path $OutDir "$Name.log"
    $mainDoc = (($OutDir -replace '\\', '/') + '/main')
    Write-Host "pdflatex $Name"
    & $Pdflatex -interaction=nonstopmode -halt-on-error -file-line-error `
        "-output-directory=$OutDir" "-jobname=$Name" `
        "\def\MainDoc{$mainDoc}\input{$Name.tex}" *> $null
    if ($LASTEXITCODE -ne 0) {
        Show-LogTail -Path $log
        if ((Test-Path -LiteralPath $log) -and
            (Select-String -LiteralPath $log -Pattern "I can't write on file" -Quiet)) {
            throw "pdflatex failed on ${Name}: output PDF is locked. Close the viewer and rerun."
        }
        throw "pdflatex failed on $Name (exit $LASTEXITCODE)"
    }
}

function Invoke-Bibtex {
    param([string]$Name)

    $arg = ($OutDir -replace '\\', '/') + "/$Name"
    $auxPath = Join-Path $OutDir "$Name.aux"
    $cited = $false
    if (Test-Path -LiteralPath $auxPath) {
        $cited = [bool](Select-String -LiteralPath $auxPath -Pattern '\\citation\{' -Quiet)
    }
    if (-not $cited) {
        Write-Host "bibtex skipped (${Name}): no citations"
        $bbl = Join-Path $OutDir "$Name.bbl"
        Set-Content -LiteralPath $bbl -Encoding ascii -Value '% no citations'
        return
    }
    Write-Host "bibtex $arg"
    & $Bibtex $arg
    $code = $LASTEXITCODE
    $blg = Join-Path $OutDir "$Name.blg"
    if ($code -ge 2) {
        Show-LogTail -Path $blg
        throw "bibtex failed on $Name (exit $code)"
    }
    if ($code -eq 1 -and (Test-Path -LiteralPath $blg)) {
        Write-Host "bibtex warnings (${Name}):"
        Select-String -LiteralPath $blg -Pattern 'Warning--|error' |
            ForEach-Object { Write-Host $_.Line }
    }
}

function Get-LogHits {
    param(
        [string]$Name,
        [string[]]$Patterns
    )
    $log = Join-Path $OutDir "$Name.log"
    return Select-String -LiteralPath $log -Pattern $Patterns
}

foreach ($name in $Docs) {
    Write-Host ""
    Write-Host "==== $name ===="
    Invoke-Pdflatex -Name $name
    Invoke-Bibtex -Name $name

    $clean = $false
    for ($pass = 1; $pass -le 4; $pass++) {
        Invoke-Pdflatex -Name $name
        $fatal = Get-LogHits -Name $name -Patterns @(
            'multiply defined',
            'multiply-defined',
            'LABELS NOT IMPORTED',
            '^! '
        )
        if ($fatal) {
            $fatal | ForEach-Object {
                Write-Host ("{0}:{1}: {2}" -f $name, $_.LineNumber, $_.Line.Trim())
            }
            throw "$name failed label import or has a LaTeX error"
        }
        $rerun = Get-LogHits -Name $name -Patterns @(
            'Rerun to get cross-references right',
            'Rerun to get it removed'
        )
        $undefined = Get-LogHits -Name $name -Patterns @(
            'undefined on input line',
            'There were undefined references'
        )
        if ((-not $rerun) -and (-not $undefined)) {
            $clean = $true
            break
        }
        Write-Host "unresolved refs, extra pass $pass ($name)"
    }
    if (-not $clean) {
        $left = Get-LogHits -Name $name -Patterns @(
            'undefined on input line',
            'There were undefined references',
            'Rerun to get cross-references right'
        )
        if ($left) {
            $left | ForEach-Object {
                Write-Host ("{0}:{1}: {2}" -f $name, $_.LineNumber, $_.Line.Trim())
            }
        }
        throw "$name still has undefined references or asks for a rerun"
    }
}

Write-Host ""
Write-Host "Done: $($Docs.Count) documents in $OutDir"
