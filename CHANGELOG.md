# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `--skip-adaptive` flag for the `qc` command to run only fixed thresholds without adaptive QC
- Params files with `adaptive: null` are automatically recognized to skip adaptive QC on `--use-params`

### Fixed
- `build_adaptive_mask` now returns an all-True mask when `interval_data` is empty, preventing all cells from being filtered out on retry failure

## [0.2.4] - 2024-08-06

### Added
- Initial public release
- Core functionality for curating single cell data
- Commands for downloading, constructing counts, QC, transformation, cell typing, and metadata harmonization
- Support for GEO dataset downloads
- Integration with scanpy for single cell analysis
- Harmony integration for batch correction
- OpenAI integration for intelligent data processing
- Docker support for containerized workflows
- Latch platform integration

### Fixed
- Package distribution now includes all Python source files and data files

### Dependencies
- scanpy 1.11.1
- scikit-image 0.25.2
- igraph 0.11.8
- harmonypy 0.0.10
- geoparse 2.0.4
- pysradb 2.2.2
- docker 7.1.0
- openai 1.81.0
- click 8.2.1
- latch >=2.63.1
- tiktoken >=0.9.0

## [0.2.3] - [Previous release date]
- Previous version history to be documented

[Unreleased]: https://github.com/yourusername/latch-curate/compare/v0.2.4...HEAD
[0.2.4]: https://github.com/yourusername/latch-curate/releases/tag/v0.2.4