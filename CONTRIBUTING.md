# Contributing to Reference Masker

We welcome contributions to Reference Masker! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your feature or fix
4. Make your changes
5. Run tests to ensure everything works
6. Submit a pull request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/reference-masker.git
cd reference-masker

# Install development dependencies
mamba create -n reference-masker-dev -c conda-forge -c bioconda snakemake biopython pyyaml matplotlib seaborn pandas numpy
conda activate reference-masker-dev

# Run tests
python test/test_workflow.py

# Validate configuration
python scripts/validate_inputs.py --config config.yaml
```

## Types of Contributions

### Bug Reports
- Use GitHub Issues to report bugs
- Include detailed description of the problem
- Provide steps to reproduce
- Include system information and log files

### Feature Requests
- Use GitHub Issues to suggest new features
- Explain the use case and benefit
- Consider implementation complexity

### Code Contributions
- Follow existing code style
- Add tests for new functionality
- Update documentation as needed
- Ensure all tests pass

## Code Style Guidelines

### Python Code
- Follow PEP 8 style guidelines
- Use descriptive variable names
- Add docstrings to functions and classes
- Include type hints where appropriate

### Snakemake Rules
- Use clear rule names
- Add comments explaining complex logic
- Include proper error handling
- Log important steps

### Shell Scripts
- Use `set -euo pipefail` for error handling
- Quote variables to prevent word splitting
- Add comments for complex operations
- Check file existence before operations

## Testing

### Running Tests
```bash
# Run all tests
python test/test_workflow.py

# Run specific test
python -m unittest test.test_workflow.TestReferenceMasker.test_validation_script
```

### Adding Tests
- Add tests for new functionality
- Test both success and failure cases
- Use descriptive test names
- Include edge cases

## Documentation

### User Documentation
- Update README.md for user-facing changes
- Include examples for new features
- Update parameter descriptions
- Add troubleshooting information

### Code Documentation
- Add docstrings to new functions
- Comment complex algorithms
- Update inline comments
- Include type hints

## Pull Request Process

1. **Create a branch**: Use descriptive branch names
   ```bash
   git checkout -b feature/add-new-validation
   git checkout -b fix/memory-leak-issue
   ```

2. **Make changes**: Follow the guidelines above

3. **Test thoroughly**:
   ```bash
   # Run tests
   python test/test_workflow.py
   
   # Test with example data
   snakemake --use-conda -j 2 --dry-run
   ```

4. **Update documentation**: Ensure README and comments are current

5. **Commit changes**: Use clear, descriptive commit messages
   ```bash
   git commit -m "Add validation for large genome files
   
   - Check file size before processing
   - Add memory estimation
   - Update user warnings"
   ```

6. **Push and create PR**:
   ```bash
   git push origin feature/add-new-validation
   ```

7. **PR Description**: Include:
   - What changes were made
   - Why the changes were necessary
   - How to test the changes
   - Any breaking changes

## Code Review Process

- All submissions require review
- Reviewers will check:
  - Code quality and style
  - Test coverage
  - Documentation updates
  - Performance implications
- Address feedback promptly
- Be open to suggestions

## Release Process

### Version Numbers
- Follow semantic versioning (MAJOR.MINOR.PATCH)
- Major: Breaking changes
- Minor: New features, backward compatible
- Patch: Bug fixes, backward compatible

### Release Checklist
- [ ] All tests pass
- [ ] Documentation updated
- [ ] Version number bumped
- [ ] Changelog updated
- [ ] Tag created
- [ ] Release notes written

## Community Guidelines

### Be Respectful
- Use welcoming and inclusive language
- Respect different viewpoints and experiences
- Accept constructive criticism gracefully
- Focus on what's best for the community

### Be Collaborative
- Help others learn and grow
- Share knowledge and resources
- Credit contributors appropriately
- Build on each other's work

## Getting Help

- **Questions**: Open a GitHub Issue with the "question" label
- **Discussions**: Use GitHub Discussions for general topics
- **Bugs**: Report via GitHub Issues with detailed information
- **Security**: Email maintainers directly for security issues

## Recognition

Contributors will be:
- Listed in the README acknowledgments
- Credited in release notes
- Invited to be maintainers for significant contributions

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to Reference Masker! 🎉