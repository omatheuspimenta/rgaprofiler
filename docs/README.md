# omatheuspimenta/rgaprofiler: Documentation

New to the pipeline? Follow the numbered steps in the main [README](../README.md#usage) —
they take you from a fresh clone to a finished test run. Then:

- [Third-party software setup](software-setup.md)
  - What license-gated software each tool needs, where to download it and where to put
    it (a required, one-time step before the first run).
- [Usage](usage.md)
  - How to run the pipeline: the samplesheet, [which setup fits your machine](usage.md#which-setup-fits-you)
    (GPU or CPU-only, small protein sets, whole proteomes, clusters), chunking, resources,
    and adapting the RGA classification to another organism.
- [Output](output.md)
  - Where each result is written (`<outdir>/<tool>/<sample>/`) and how to read it.

For maintainers:

- [Publishing the Docker images](publishing-docker-images.md) — building, tagging and
  publishing the pipeline's container images, including how to release an updated image.
- [Contributing](CONTRIBUTING.md)
