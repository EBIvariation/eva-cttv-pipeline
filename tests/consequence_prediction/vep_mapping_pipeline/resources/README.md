# Tests

The pipeline includes a simple integration test. It runs a sample VCF of 2,000 records through VEP, using a wrapper script, and checks that the run has been successful and that the output records match the expected ones. Further details are explained in [`.travis.yml`](/.travis.yml).