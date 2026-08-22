FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages: OpenJDK 11, Perl, Python3, procps, and libpcre3
RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    perl \
    python3 \
    procps \
    libpcre3 \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME="/usr/lib/jvm/java-11-openjdk-amd64"
ENV PATH="$JAVA_HOME/bin:${PATH}"

# We don't need ENTRYPOINT and we don't download the database here!