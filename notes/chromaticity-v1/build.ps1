# Build chromaticity note: pdflatex + bibtex (TeX Live 2024), output build\chromaticity.pdf
$ErrorActionPreference = 'Stop'
$TexLive = 'C:\texlive\2024\bin\windows'
$env:PATH = "$TexLive;$env:PATH"; $env:BIBINPUTS = '.;'
Set-Location -LiteralPath $PSScriptRoot
New-Item -ItemType Directory -Force -Path build | Out-Null
function Tex { & pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build -jobname=chromaticity main.tex *> $null; if ($LASTEXITCODE -ne 0) { Get-Content build\chromaticity.log -Tail 40; throw 'pdflatex failed' } }
Tex; & bibtex build/chromaticity; Tex; Tex; Tex
$bad = Select-String build\chromaticity.log -Pattern 'undefined on input line','There were undefined references','Rerun to get','^! '
if ($bad) { $bad | % Line; throw 'unresolved references' }
Write-Host 'Done: build\chromaticity.pdf'
