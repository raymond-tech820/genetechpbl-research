param(
    [switch]$CpuOnly
)

$ErrorActionPreference = "Stop"
$ReleaseRoot = Split-Path -Parent (Split-Path $MyInvocation.MyCommand.Path)
$PipelineRoot = Join-Path $ReleaseRoot "research_pipeline"
$VenvPath = Join-Path $ReleaseRoot ".venv"
$PythonCommand = "python"

if (-not (Test-Path -LiteralPath $VenvPath)) {
    & $PythonCommand -m venv $VenvPath
}

$VenvPython = Join-Path $VenvPath "Scripts\\python.exe"
& $VenvPython -m pip install --upgrade pip setuptools wheel

if ($CpuOnly) {
    & $VenvPython -m pip install "torch==2.13.0" --index-url https://download.pytorch.org/whl/cpu
} else {
    & $VenvPython -m pip install "torch==2.13.0" --index-url https://download.pytorch.org/whl/cu130
}

& $VenvPython -m pip install -r (Join-Path $PipelineRoot "requirements.lock")
& $VenvPython -m pip install --no-deps -e $PipelineRoot
& $VenvPython -m pytest (Join-Path $PipelineRoot "tests")
