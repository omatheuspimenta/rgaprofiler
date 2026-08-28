FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages: OpenJDK 11, Perl, Python3, procps, and libpcre3.
# libgomp1 (GNU OpenMP runtime) is required by some member-database binaries
# (e.g. bin/cdd/rpsblast) even though InterProScan's own docs don't call it
# out explicitly — without it those analyses fail with "libgomp.so.1: cannot
# open shared object file".
RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    perl \
    python3 \
    procps \
    libpcre3 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME="/usr/lib/jvm/java-11-openjdk-amd64"
ENV PATH="$JAVA_HOME/bin:${PATH}"

# We don't need ENTRYPOINT and we don't download the database here!

# Run as a non-root user by default. Nextflow's docker profile overrides this via
# `-u $(id -u):$(id -g)` at task execution time regardless (needed anyway since
# InterProScan writes temp/work dirs inside the user-mounted --interproscan_db
# directory), but this keeps the image itself sane for anyone running it directly.
RUN useradd --create-home --shell /bin/bash interproscan
USER interproscan
WORKDIR /data