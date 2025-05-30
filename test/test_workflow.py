#!/usr/bin/env python3
"""
Test suite for Reference Masker workflow
"""

import os
import sys
import unittest
import tempfile
import shutil
from pathlib import Path
import subprocess

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

class TestReferenceMasker(unittest.TestCase):
    """Test cases for Reference Masker workflow"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        
        # Create test viral database
        self.viral_db = os.path.join(self.test_dir, "test_viruses.fasta")
        with open(self.viral_db, 'w') as f:
            f.write(">test_virus_1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
            f.write(">test_virus_2\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n")
        
        # Create test reference genome
        self.ref_genome = os.path.join(self.test_dir, "test_ref.fasta")
        with open(self.ref_genome, 'w') as f:
            f.write(">test_chromosome\n")
            # Include some sequence that matches viral database
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
            # Add some random sequence
            f.write("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            f.write("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")
            f.write("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
        
        # Create test config
        self.config_file = os.path.join(self.test_dir, "test_config.yaml")
        with open(self.config_file, 'w') as f:
            f.write(f"""
viral_database: "{self.viral_db}"
reference_genomes:
  - "{self.ref_genome}"

output_dir: "{os.path.join(self.test_dir, 'results')}"
temp_dir: "{os.path.join(self.test_dir, 'temp')}"
log_dir: "{os.path.join(self.test_dir, 'logs')}"

shred:
  overlap: 20
  median_length: 40
  variance: 5
  min_length: 30

mapping:
  min_identity: 90
  max_indel: 2

masking:
  entropy_threshold: 0.5

threads: 2
""")
    
    def tearDown(self):
        """Clean up test fixtures"""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)
    
    def test_validation_script(self):
        """Test input validation script"""
        # Import validation functions
        from validate_inputs import validate_config_file, validate_fasta_file
        
        # Test config validation
        errors, warnings = validate_config_file(self.config_file)
        self.assertEqual(len(errors), 0, f"Config validation failed: {errors}")
        
        # Test FASTA validation
        errors, warnings = validate_fasta_file(self.viral_db, "Test viral database")
        self.assertEqual(len(errors), 0, f"Viral database validation failed: {errors}")
        
        errors, warnings = validate_fasta_file(self.ref_genome, "Test reference")
        self.assertEqual(len(errors), 0, f"Reference genome validation failed: {errors}")
    
    def test_file_existence(self):
        """Test that required files exist"""
        self.assertTrue(os.path.exists(self.viral_db))
        self.assertTrue(os.path.exists(self.ref_genome))
        self.assertTrue(os.path.exists(self.config_file))
    
    def test_fasta_format(self):
        """Test FASTA file format"""
        # Test viral database
        with open(self.viral_db, 'r') as f:
            content = f.read()
            self.assertTrue(content.startswith('>'))
            self.assertIn('ATCG', content)
        
        # Test reference genome
        with open(self.ref_genome, 'r') as f:
            content = f.read()
            self.assertTrue(content.startswith('>'))
            self.assertIn('ATCG', content)
    
    def test_config_format(self):
        """Test configuration file format"""
        import yaml
        
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Test required sections
        required_sections = ['viral_database', 'reference_genomes', 'output_dir', 'temp_dir', 'log_dir']
        for section in required_sections:
            self.assertIn(section, config)
        
        # Test parameter sections
        self.assertIn('shred', config)
        self.assertIn('mapping', config)
        self.assertIn('masking', config)
    
    def test_snakemake_dry_run(self):
        """Test Snakemake dry run"""
        # Change to test directory
        os.chdir(self.test_dir)
        
        # Copy necessary files to test directory
        workflow_dir = os.path.dirname(os.path.dirname(__file__))
        
        # Copy Snakefile
        shutil.copy(os.path.join(workflow_dir, 'Snakefile'), self.test_dir)
        
        # Copy environment files
        env_dir = os.path.join(self.test_dir, 'envs')
        os.makedirs(env_dir, exist_ok=True)
        shutil.copy(os.path.join(workflow_dir, 'envs', 'bbmap.yaml'), env_dir)
        shutil.copy(os.path.join(workflow_dir, 'envs', 'plotting.yaml'), env_dir)
        
        # Copy scripts
        scripts_dir = os.path.join(self.test_dir, 'scripts')
        os.makedirs(scripts_dir, exist_ok=True)
        shutil.copy(os.path.join(workflow_dir, 'scripts', 'generate_report.py'), scripts_dir)
        
        try:
            # Run snakemake dry run
            result = subprocess.run([
                'snakemake', '--dry-run', '--configfile', 'test_config.yaml'
            ], capture_output=True, text=True, timeout=30)
            
            # Check that dry run succeeded
            self.assertEqual(result.returncode, 0, 
                           f"Snakemake dry run failed: {result.stderr}")
            
        except subprocess.TimeoutExpired:
            self.fail("Snakemake dry run timed out")
        except FileNotFoundError:
            self.skipTest("Snakemake not available in test environment")

class TestValidationScript(unittest.TestCase):
    """Test cases for validation script specifically"""
    
    def test_invalid_config(self):
        """Test validation with invalid configuration"""
        from validate_inputs import validate_config_file
        
        # Test missing file
        errors, warnings = validate_config_file("nonexistent.yaml")
        self.assertGreater(len(errors), 0)
    
    def test_invalid_fasta(self):
        """Test validation with invalid FASTA file"""
        from validate_inputs import validate_fasta_file
        
        # Test missing file
        errors, warnings = validate_fasta_file("nonexistent.fasta")
        self.assertGreater(len(errors), 0)

def run_tests():
    """Run all tests"""
    # Set up test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test cases
    suite.addTests(loader.loadTestsFromTestCase(TestReferenceMasker))
    suite.addTests(loader.loadTestsFromTestCase(TestValidationScript))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)