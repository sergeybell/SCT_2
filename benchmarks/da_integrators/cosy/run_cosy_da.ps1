# run_cosy_da.ps1 - generate da_ops_o<N>.fox from the template and run COSY for N = 3, 5, 7, 10
$cosy = 'C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe'
Set-Location $PSScriptRoot; New-Item -ItemType Directory output -Force | Out-Null
foreach ($n in 3, 5, 7, 10) {
  $s = (Get-Content da_ops_template.fox -Raw).Replace('@NO@', "$n")
  [IO.File]::WriteAllText("$PWD\da_ops_o$n.fox", $s)
  & $cosy "da_ops_o$n.fox" > "output\da_ops_o$n.log" 2>&1
  Write-Host "order $n done"; Get-Content "output\time_o$n.txt"
}