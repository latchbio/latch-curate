build:
  rm -rf dist
  uv build

local-install:
	uv pip uninstall .
	uv pip install -e .
