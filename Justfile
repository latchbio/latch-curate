build:
  rm -rf dist
  uv build

local-install:
  uv pip uninstall .
  uv pip install -e .

publish:
  uv publish --token $(<credentials/pypi_token)
  rm -rf dist
