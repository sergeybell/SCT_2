# run_cosy_ode.ps1 - Henon-Heiles (FOX RK4) and lone elements MQ / MS for N = 3, 5, 7, 10
$cosy = 'C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe'
Set-Location $PSScriptRoot; New-Item -ItemType Directory output -Force | Out-Null
foreach ($n in 3, 5, 7, 10) { foreach ($t in 'hh', 'elem') {
  $s = (Get-Content "${t}_template.fox" -Raw).Replace('@NO@', "$n")
  [IO.File]::WriteAllText("$PWD\${t}_o$n.fox", $s)
  & $cosy "${t}_o$n.fox" > "output\${t}_o$n.log" 2>&1
  Write-Host "$t order $n done" } }