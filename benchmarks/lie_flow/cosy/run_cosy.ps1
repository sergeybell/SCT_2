# run_cosy.ps1 - lone DL / MQ / MS / MH maps (all 6 rows) for N = 3, 5, 7, 10
$cosy = 'C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe'
Set-Location $PSScriptRoot; New-Item -ItemType Directory output -Force | Out-Null
foreach ($n in 3, 5, 7, 10) {
  $s = (Get-Content elem6_template.fox -Raw).Replace('@NO@', "$n")
  [IO.File]::WriteAllText("$PWD\elem6_o$n.fox", $s)
  & $cosy "elem6_o$n.fox" > "output\elem6_o$n.log" 2>&1
  Write-Host "order $n done"; Get-Content "output\elem_time_o$n.txt" }