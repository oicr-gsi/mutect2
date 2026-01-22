# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-01-23
- [GRD-1011](https://jira.oicr.on.ca/browse/GRD-1011) - Added new wdl task to calculate cache metrics in the workflow.

## [1.1.3] - 2025-11-28
- identical to 1.1.1, somehow 1.1.1 has already been built without any record or being visible in jenkins.
- Incrementing version for deployment

## [1.1.1] - 2025-11-21
### Added
- Baked-in gatk module set to "gatk/4.2.6.1"

## [Unreleased] - 2025-04-23
- [GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - Expanded built-in documentation (metadata changes only).

## [1.1.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to metadata only)

## [1.0.9] - 2024-02-12
### Added
- Added minMemory parameter for runMutect2 task to better handle RAM allocation for small chromosomes

## [1.0.8] - 2024-01-20
### Added
- GRD-728 Introducing RAM scaling by chromosome

## [Unreleased] - 2021-11-25
### Fixed
- GP-2883 Making RT tests more robust

## [1.0.7] - 2023-10-12
### Added
- GCGI-995 Add Gnomad AF to hg38 input

## [1.0.6] - 2023-07-26
### Added
- GRD-674. Adding normal_only alias as one step creating PON  

## [1.0.5] - 2023-05-29
### Added
- Grd 492 improvements.

## [1.0.4] - 2021-08-03
### Added
- Use String path for interval bed file.

## [1.0.3] - 2020-06-22
### Changed
- Rename from mutect2GATK4 to mutect2

## [Unreleased]
### Fixed
- Hotfix to expose -pon and --germline-resource

## [1.0.0] - 2021-05-31
### Added
- Add correct workflow names to vidarrbuild.
