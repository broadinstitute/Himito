# GitHub Actions Workflows

This directory contains GitHub Actions workflows for automated testing and building of Himito.

## Workflows

### `ci.yml` - Continuous Integration
- **Triggers**: Push/PR to main/master/develop branches
- **Tests**: 
  - Code formatting (cargo fmt)
  - Linting (cargo clippy)
  - Unit tests
  - Build verification
  - Command-line interface tests
- **Platforms**: Ubuntu and macOS
- **Rust versions**: stable and beta

### `test-small.yml` - Quick Integration Tests
- **Triggers**: Push/PR to main/master branches
- **Tests**:
  - Command-line help verification
  - Integration tests with minimal synthetic data
  - Evaluate command functionality
- **Purpose**: Fast feedback on basic functionality

### `release.yml` - Release Builds
- **Triggers**: GitHub releases or manual dispatch
- **Builds**: Release binaries for Linux, macOS, and Windows
- **Outputs**: Platform-specific archives (.tar.gz, .zip)

## Running Tests Locally

To run tests locally before pushing:

```bash
# Format code
cargo fmt --all

# Build release
cargo build --release

# Test help commands
./target/release/haplograph --help
./target/release/haplograph haplograph --help
./target/release/haplograph dev-tools --help
```

## Adding New Tests

To add new integration tests:

1. Add test data to `test/` directory (if needed)
2. Add test steps to `test-small.yml` or `ci.yml`
3. Ensure tests are idempotent and don't require external services
4. Use conditional checks for optional dependencies

## Notes

- Tests that require large BAM files or external data are skipped if files are not available
- All workflows cache Rust dependencies for faster builds
- Release builds are only created on actual releases, not on every push

